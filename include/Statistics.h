#pragma once
#include <vector>

namespace DistTool {

/**
 * Statistical utilities for SOAP descriptor analysis.
 *
 * Implements the metrics described in the paper:
 *   - Mean descriptor vector  q̄(T)  from a thermalized reference sample.
 *   - Euclidean distance  d^i = ||q̃^i − q̄(T)||  (Eq. 5).
 *   - Chi-distribution probability model  P(q̃^i | T)  (Eq. 6).
 */
class Statistics {
public:
    // Arithmetic mean of a set of vectors.
    static std::vector<double> mean(
        const std::vector<std::vector<double>>& data);

    // Euclidean distance  ||a − b||.
    static double euclidean(const std::vector<double>& a,
                             const std::vector<double>& b);

    // Chi-distribution probability model (FaVaD Eq. 6):
    //   P(d | k, σ) ∝ d^(k−2) · exp(−d²/(2σ²))
    // Normalised to 1 at the distribution mode d_peak = σ·√(k−2)  (k > 2).
    // k is the effective number of active DV components; σ is a scale fitted
    // from the reference distance distribution.
    static double chiProbability(double d, double k, double sigma);

    // Fit chi-distribution parameters (k, sigma) from a set of distances
    // using method of moments:
    //   k     = 2·<d²>² / Var(d²)
    //   sigma = sqrt(<d²> / k)
    static void fitChiParams(const std::vector<double>& dists,
                             double& k_out, double& sigma_out);

    // Non-parametric anomaly score: the empirical CDF of `d` within a sorted
    // (ascending) reference distance sample — i.e. the fraction of reference
    // atoms whose distance is ≤ d. Returns a value in [0, 1] that behaves as a
    // percentile rank: ~0 for a typical environment, →1 for an outlier. This is
    // robust to non-chi-shaped distributions (e.g. chemically disordered alloys
    // where the chi method-of-moments fit degenerates). Empty sample → 0.
    static double empiricalCdf(const std::vector<double>& sorted_asc, double d);

    // The p-th percentile of a sorted (ascending) sample, p ∈ [0,1], with
    // linear interpolation between order statistics. Empty sample → 0.
    static double percentile(const std::vector<double>& sorted_asc, double p);
};

} // namespace DistTool
