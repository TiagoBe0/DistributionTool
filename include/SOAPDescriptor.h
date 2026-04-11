#pragma once
#include "AtomData.h"
#include "RadialBasis.h"
#include "SphericalHarmonics.h"
#include <array>
#include <vector>

namespace DistTool {

/**
 * SOAP (Smooth Overlap of Atomic Positions) descriptor.
 *
 * For each atom i the descriptor vector (DV) is the power spectrum:
 *
 *   c_{nlm}^i = Σ_j  φ_n(r_{ij}) · Y_l^m( r̂_{ij} )
 *
 *   p_{nn'l}^i = π · sqrt(8/(2l+1)) · Σ_{m=−l}^{l}  c_{nlm}^i · c_{n'lm}^i
 *
 * Only the upper-triangular (n ≤ n') terms are stored, giving
 *   DV size = (n_max·(n_max+1)/2) · (l_max+1)
 *
 * The final DV is normalised to unit length: q̃^i = p^i / |p^i|.
 *
 * Reference: Bartók et al., Phys. Rev. B 87 (2013) 184115.
 */
struct SOAPParams {
    int    n_max     = 9;      // Radial basis functions
    int    l_max     = 9;      // Max angular momentum
    double r_cut     = 5.0;   // Cutoff radius [Å]
    double sigma     = -1.0;  // Gaussian width (-1 = auto = spacing/1)
    bool   normalize = true;  // Normalise DV to unit length

    // Number of DV components
    int dvSize() const {
        return (n_max * (n_max + 1) / 2) * (l_max + 1);
    }
};

class SOAPDescriptor {
public:
    explicit SOAPDescriptor(const SOAPParams& p);

    // Compute the DV for a single atom given its neighbour displacement vectors.
    // neighbors: list of (Δx, Δy, Δz) already in Å with PBC applied.
    std::vector<double> compute(
        const std::vector<std::array<double,3>>& neighbors) const;

    // Compute DVs for every atom in the frame (applies PBC internally).
    void computeAll(Frame& frame) const;

    const SOAPParams& params() const { return params_; }
    int dvSize() const { return params_.dvSize(); }

private:
    SOAPParams         params_;
    SphericalHarmonics sph_;
    RadialBasis        rad_;
    std::vector<double> ps_factors_;  // π·√(8/(2l+1)) per l — pre-computed

    // Inner computation kernel — caller owns all buffers (no heap allocs inside).
    void computeWithBuffers(
        const std::vector<std::array<double,3>>& neighbors,
        double* phi_buf,    // [n_max]
        double* ylm_buf,    // [(l_max+1)²]
        double* c_buf,      // [n_max × (l_max+1)²], zeroed before call
        std::vector<double>& dv_out) const;

    // Fallback O(N) brute-force neighbour list (used if cell list unavailable).
    std::vector<std::array<double,3>> neighborList(
        const Frame& frame, int i) const;

    static std::array<double,3> minImage(
        std::array<double,3> dr, const SimBox& box);
};

} // namespace DistTool
