#pragma once
#include <vector>

namespace DistTool {

/**
 * Gaussian radial basis functions used in the SOAP descriptor.
 *
 * n_max equidistant Gaussian centres along (0, r_cut]:
 *   r_n = (n+1) · r_cut / (n_max+1)
 *
 * φ_n(r) = exp( −(r − r_n)² / (2σ²) ) · f_cut(r)
 *
 * f_cut(r) = 0.5·(1 + cos(π·r/r_cut))  for r < r_cut, else 0
 *
 * Performance: all sigma/spacing constants pre-computed at construction.
 * Use computeInto() to write into a caller-owned buffer — no heap allocation.
 */
class RadialBasis {
public:
    RadialBasis(int n_max, double r_cut, double sigma = -1.0);

    // Write all φ_n(r) into pre-allocated array out[n_max].
    // Returns false (and zeros out[]) if r >= r_cut.  HOT PATH — no malloc.
    bool computeInto(double r, double* out) const;

    // Convenience wrapper returning a new vector.
    std::vector<double> compute(double r) const;

    double cutoff(double r) const;

    int    nMax()  const { return n_max_; }
    double rCut()  const { return r_cut_; }
    double sigma() const { return sigma_; }

private:
    int    n_max_;
    double r_cut_;
    double sigma_;
    double inv2s2_;           // 1 / (2σ²)  — pre-computed
    std::vector<double> centres_;
};

} // namespace DistTool
