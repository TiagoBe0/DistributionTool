#include "SOAPDescriptor.h"
#include "CellList.h"
#include <cmath>
#include <numeric>
#include <stdexcept>

#ifdef USE_OPENMP
#  include <omp.h>
#endif

namespace DistTool {

SOAPDescriptor::SOAPDescriptor(const SOAPParams& p)
    : params_(p)
    , sph_(p.l_max)
    , rad_(p.n_max, p.r_cut, p.sigma)
{
    // Pre-compute power-spectrum prefactors  π·√(8/(2l+1))  — one per l value.
    // Avoids recomputing std::sqrt inside the tight inner loop for every atom.
    ps_factors_.resize(p.l_max + 1);
    for (int l = 0; l <= p.l_max; ++l)
        ps_factors_[l] = M_PI * std::sqrt(8.0 / (2.0 * l + 1.0));
}

// ── Minimum-image convention ──────────────────────────────────────────────────
std::array<double,3> SOAPDescriptor::minImage(
    std::array<double,3> dr, const SimBox& box)
{
    const double L[3] = { box.lx(), box.ly(), box.lz() };
    for (int d = 0; d < 3; ++d)
        if (box.periodic[d] && L[d] > 0.0)
            dr[d] -= L[d] * std::round(dr[d] / L[d]);
    return dr;
}

// ── Fallback brute-force neighbour list (used only when cell list unavailable)
std::vector<std::array<double,3>> SOAPDescriptor::neighborList(
    const Frame& frame, int i) const
{
    const Atom& ai = frame.atoms[i];
    const double rc2 = params_.r_cut * params_.r_cut;
    std::vector<std::array<double,3>> nb;
    nb.reserve(64);
    for (int j = 0; j < frame.size(); ++j) {
        if (j == i) continue;
        const Atom& aj = frame.atoms[j];
        std::array<double,3> dr = {aj.x - ai.x, aj.y - ai.y, aj.z - ai.z};
        dr = minImage(dr, frame.box);
        if (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] < rc2)
            nb.push_back(dr);
    }
    return nb;
}

// ── Core SOAP computation — HOT PATH ─────────────────────────────────────────
// All work buffers are passed in by the caller so no heap allocation happens
// inside this function.  `phi_buf` size = n_max, `ylm_buf` size = nYlm,
// `c_buf` size = n_max * nYlm (zeroed by caller before first use).
void SOAPDescriptor::computeWithBuffers(
    const std::vector<std::array<double,3>>& neighbors,
    double* phi_buf,    // [n_max]
    double* ylm_buf,    // [nYlm]
    double* c_buf,      // [n_max * nYlm], must be zeroed before call
    std::vector<double>& dv_out) const
{
    const int nmax = params_.n_max;
    const int lmax = params_.l_max;
    const int nYlm = (lmax + 1) * (lmax + 1);

    for (const auto& dr : neighbors) {
        const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        const double r  = std::sqrt(r2);
        if (r < 1e-12 || r >= params_.r_cut) continue;

        rad_.computeInto(r, phi_buf);
        sph_.computeInto(dr[0], dr[1], dr[2], ylm_buf);

        // Accumulate c[n][lm] += phi_n * ylm  (flat layout: c_buf[n*nYlm + lm])
        for (int n = 0; n < nmax; ++n) {
            const double phi_n = phi_buf[n];
            double* cn = c_buf + n * nYlm;
            for (int lm = 0; lm < nYlm; ++lm)
                cn[lm] += phi_n * ylm_buf[lm];
        }
    }

    // Power spectrum  p_{nn'l} = factor_l · Σ_m c_{nlm}·c_{n'lm}
    // Only upper triangle  n ≤ n'
    const int dv_size = params_.dvSize();
    dv_out.resize(dv_size);

    int idx = 0;
    for (int n = 0; n < nmax; ++n) {
        const double* cn = c_buf + n * nYlm;
        for (int np = n; np < nmax; ++np) {
            const double* cnp = c_buf + np * nYlm;
            for (int l = 0; l <= lmax; ++l) {
                double sum = 0.0;
                const int base = l * l + l;
                for (int m = -l; m <= l; ++m)
                    sum += cn[base + m] * cnp[base + m];
                dv_out[idx++] = ps_factors_[l] * sum;
            }
        }
    }

    // Normalise to unit length  q̃ = p / |p|
    if (params_.normalize) {
        double norm2 = 0.0;
        for (double v : dv_out) norm2 += v * v;
        if (norm2 > 1e-20) {
            const double inv_norm = 1.0 / std::sqrt(norm2);
            for (auto& v : dv_out) v *= inv_norm;
        }
    }
}

// ── Public single-atom interface (allocates its own buffers) ──────────────────
std::vector<double> SOAPDescriptor::compute(
    const std::vector<std::array<double,3>>& neighbors) const
{
    const int nYlm = (params_.l_max + 1) * (params_.l_max + 1);
    std::vector<double> phi_buf(params_.n_max);
    std::vector<double> ylm_buf(nYlm);
    std::vector<double> c_buf(params_.n_max * nYlm, 0.0);
    std::vector<double> dv;
    computeWithBuffers(neighbors, phi_buf.data(), ylm_buf.data(), c_buf.data(), dv);
    return dv;
}

// ── Compute DVs for all atoms: cell list + OpenMP + per-thread buffers ────────
void SOAPDescriptor::computeAll(Frame& frame) const {
    const int N    = frame.size();
    const int nmax = params_.n_max;
    const int nYlm = (params_.l_max + 1) * (params_.l_max + 1);

    // Build the cell list once — O(N), avoids O(N²) brute-force neighbour search
    CellList cl;
    cl.build(frame, params_.r_cut);

#ifdef USE_OPENMP
    #pragma omp parallel
    {
        // Per-thread work buffers — allocated ONCE per thread, not per atom.
        // This eliminates millions of malloc/free calls from the inner loop.
        std::vector<double> phi_buf(nmax);
        std::vector<double> ylm_buf(nYlm);
        std::vector<double> c_buf(nmax * nYlm);
        std::vector<double> dv;

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            auto nb = cl.neighbours(i);

            // Zero the coefficient buffer before accumulation
            std::fill(c_buf.begin(), c_buf.end(), 0.0);

            computeWithBuffers(nb, phi_buf.data(), ylm_buf.data(),
                               c_buf.data(), dv);
            frame.atoms[i].dv = dv;
        }
    }
#else
    std::vector<double> phi_buf(nmax);
    std::vector<double> ylm_buf(nYlm);
    std::vector<double> c_buf(nmax * nYlm);
    std::vector<double> dv;

    for (int i = 0; i < N; ++i) {
        auto nb = cl.neighbours(i);
        std::fill(c_buf.begin(), c_buf.end(), 0.0);
        computeWithBuffers(nb, phi_buf.data(), ylm_buf.data(), c_buf.data(), dv);
        frame.atoms[i].dv = dv;
    }
#endif
}

} // namespace DistTool
