#include "SphericalHarmonics.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DistTool {

SphericalHarmonics::SphericalHarmonics(int l_max) : l_max_(l_max) {
    if (l_max < 0) throw std::invalid_argument("l_max must be >= 0");

    // Pre-compute normalization table once (no factorial loops at runtime)
    const int nc = (l_max + 1) * (l_max + 1);
    norm_.resize(nc);
    for (int l = 0; l <= l_max; ++l) {
        norm_[idx(l, 0)] = Nlm(l, 0);                          // m = 0
        for (int m = 1; m <= l; ++m) {
            double n = std::sqrt(2.0) * Nlm(l, m);
            norm_[idx(l,  m)] = n;   // used for cos part (Y_l^{+m})
            norm_[idx(l, -m)] = n;   // used for sin part (Y_l^{-m})
        }
    }
}

// ── N_{l,m} = sqrt( (2l+1)/(4π) · (l−m)!/(l+m)! ) ──────────────────────────
double SphericalHarmonics::Nlm(int l, int m) {
    double val = (2.0 * l + 1.0) / (4.0 * M_PI);
    for (int i = l - m + 1; i <= l + m; ++i) val /= static_cast<double>(i);
    return std::sqrt(val);
}

// ── Hot-path: write all Y_l^m into pre-allocated array ───────────────────────
//
// Optimisations over a naïve implementation:
//
//  1. Single-pass Plm recurrence (Bonnet):
//       P_{m+1}^{m+1} = -(2m+1)·sinθ · P_m^m          (diagonal step)
//       P_{m+1}^m     = x·(2m+1) · P_m^m               (superdiag step)
//       P_l^m = [x·(2l-1)·P_{l-1}^m - (l+m-1)·P_{l-2}^m] / (l-m)
//     Each P_l^m is computed exactly once: O(l_max²) instead of O(l_max³).
//
//  2. Trig recurrence replaces atan2 + cos/sin(mφ) calls:
//       cos((m+1)φ) = cos(mφ)·cos(φ) - sin(mφ)·sin(φ)
//     cos(φ) and sin(φ) are derived from dx/rxy, dy/rxy — no atan2.
//
void SphericalHarmonics::computeInto(double dx, double dy, double dz, double* out) const {
    const int nYlm = (l_max_ + 1) * (l_max_ + 1);
    std::fill(out, out + nYlm, 0.0);

    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (r < 1e-12) {
        out[0] = norm_[0];   // Y_0^0 = N_{0,0} · P_0^0
        return;
    }

    const double x   = dz / r;                    // cos θ
    const double rxy = std::sqrt(dx*dx + dy*dy);  // r·sin θ
    const double st  = rxy / r;                   // sin θ ≥ 0

    // cos φ, sin φ derived from Cartesian coords — avoids atan2 entirely.
    // When sin θ ≈ 0 (poles), all Y_l^m with m≠0 vanish anyway (P_l^m(±1)=0),
    // so the fallback values (1, 0) are never multiplied by anything nonzero.
    const double cp = (rxy > 1e-12) ? dx / rxy : 1.0;
    const double sp = (rxy > 1e-12) ? dy / rxy : 0.0;

    // ── m = 0 : Y_l^0 = N_{l,0} · P_l^0(cosθ) ──────────────────────────────
    {
        double pm_prev = 1.0;  // P_0^0
        double pm_cur  = x;    // P_1^0
        out[idx(0, 0)] = norm_[idx(0, 0)] * pm_prev;
        if (l_max_ >= 1)
            out[idx(1, 0)] = norm_[idx(1, 0)] * pm_cur;
        for (int l = 2; l <= l_max_; ++l) {
            const double pm_next = (x*(2.0*l - 1.0)*pm_cur - (l - 1.0)*pm_prev) / l;
            pm_prev = pm_cur;
            pm_cur  = pm_next;
            out[idx(l, 0)] = norm_[idx(l, 0)] * pm_cur;
        }
    }

    // ── m ≥ 1 : Y_l^{±m} = N_{l,m}·P_l^m(cosθ)·{cos,sin}(mφ) ─────────────
    // pmm  = P_m^m,    initialised to P_1^1 = -sinθ.
    // trig = (cos(mφ), sin(mφ)), initialised to m = 1.
    double pmm      = -st;
    double cos_mphi = cp;
    double sin_mphi = sp;

    for (int m = 1; m <= l_max_; ++m) {
        // Store Y_m^m
        {
            const double NP = norm_[idx(m, m)] * pmm;
            out[idx(m,  m)] = NP * cos_mphi;
            out[idx(m, -m)] = NP * sin_mphi;
        }

        // Advance along l with the Bonnet two-term recurrence
        double pm_prev = pmm;
        double pm_cur  = x * (2.0*m + 1.0) * pmm;   // P_{m+1}^m

        if (m + 1 <= l_max_) {
            const double NP = norm_[idx(m+1, m)] * pm_cur;
            out[idx(m+1,  m)] = NP * cos_mphi;
            out[idx(m+1, -m)] = NP * sin_mphi;
        }

        for (int l = m + 2; l <= l_max_; ++l) {
            const double pm_next =
                (x*(2.0*l - 1.0)*pm_cur - (l + m - 1.0)*pm_prev) / (l - m);
            pm_prev = pm_cur;
            pm_cur  = pm_next;
            const double NP = norm_[idx(l, m)] * pm_cur;
            out[idx(l,  m)] = NP * cos_mphi;
            out[idx(l, -m)] = NP * sin_mphi;
        }

        // Advance trig recurrence for m+1
        const double cos_next = cos_mphi * cp - sin_mphi * sp;
        const double sin_next = sin_mphi * cp + cos_mphi * sp;
        cos_mphi = cos_next;
        sin_mphi = sin_next;

        // Advance diagonal: P_{m+1}^{m+1} = -(2m+1)·sinθ·P_m^m
        pmm *= -(2.0*m + 1.0) * st;
    }
}

// ── Convenience wrapper ───────────────────────────────────────────────────────
std::vector<double> SphericalHarmonics::compute(double dx, double dy, double dz) const {
    std::vector<double> out(nComponents());
    computeInto(dx, dy, dz, out.data());
    return out;
}

} // namespace DistTool
