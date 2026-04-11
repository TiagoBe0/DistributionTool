#include "RadialBasis.h"
#include <cmath>
#include <stdexcept>

namespace DistTool {

RadialBasis::RadialBasis(int n_max, double r_cut, double sigma)
    : n_max_(n_max), r_cut_(r_cut)
{
    if (n_max <= 0) throw std::invalid_argument("n_max must be > 0");
    if (r_cut <= 0) throw std::invalid_argument("r_cut must be > 0");

    const double spacing = r_cut / (n_max + 1.0);
    centres_.resize(n_max);
    for (int n = 0; n < n_max; ++n)
        centres_[n] = (n + 1.0) * spacing;

    sigma_   = (sigma > 0.0) ? sigma : spacing;
    inv2s2_  = 1.0 / (2.0 * sigma_ * sigma_);   // pre-computed once
}

double RadialBasis::cutoff(double r) const {
    if (r <= 0.0)    return 1.0;
    if (r >= r_cut_) return 0.0;
    return 0.5 * (1.0 + std::cos(M_PI * r / r_cut_));
}

// ── Hot-path — no heap allocation ────────────────────────────────────────────
bool RadialBasis::computeInto(double r, double* out) const {
    if (r >= r_cut_) {
        for (int n = 0; n < n_max_; ++n) out[n] = 0.0;
        return false;
    }
    const double fc = cutoff(r);
    for (int n = 0; n < n_max_; ++n) {
        const double dr = r - centres_[n];
        out[n] = fc * std::exp(-dr * dr * inv2s2_);
    }
    return true;
}

std::vector<double> RadialBasis::compute(double r) const {
    std::vector<double> phi(n_max_);
    computeInto(r, phi.data());
    return phi;
}

} // namespace DistTool
