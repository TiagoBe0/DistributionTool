#pragma once
#include <cmath>
#include <vector>

namespace DistTool {

/**
 * Real (tesseral) spherical harmonics  Y_l^m(θ, φ).
 *
 * Convention (same as QUIP/libatoms):
 *   m = 0  :  Y_l^0  = N_{l,0} · P_l^0(cos θ)
 *   m > 0  :  Y_l^m  = √2 · N_{l,m} · P_l^m(cos θ) · cos(m φ)
 *   m < 0  :  Y_l^m  = √2 · N_{l,|m|} · P_l^|m|(cos θ) · sin(|m| φ)
 *
 * Storage order (index = l² + l + m):
 *   (0,0), (1,−1),(1,0),(1,1), (2,−2),...
 *
 * Performance: N_{l,m} normalization constants are pre-computed at construction
 * and stored in a flat table — no per-call factorial loops.
 */
class SphericalHarmonics {
public:
    explicit SphericalHarmonics(int l_max);

    // Write all Y_l^m values into pre-allocated array out[nComponents()].
    // Pass the Cartesian direction directly; normalisation is internal.
    // HOT PATH — no heap allocations.
    void computeInto(double dx, double dy, double dz, double* out) const;

    // Convenience wrapper returning a new vector.
    std::vector<double> compute(double dx, double dy, double dz) const;

    int lMax()        const { return l_max_; }
    int nComponents() const { return (l_max_ + 1) * (l_max_ + 1); }

    static int idx(int l, int m) { return l * l + l + m; }

private:
    int l_max_;

    // Pre-computed table: norm_[idx(l,m)] = N_{l,m}       for m = 0
    //                                      = √2 · N_{l,m}  for m > 0
    std::vector<double> norm_;

    // Normalization  N_{l,m} = sqrt( (2l+1)/(4π) · (l−m)!/(l+m)! )
    static double Nlm(int l, int m);
};

} // namespace DistTool
