#pragma once
#include <Eigen/Dense>
#include <vector>

namespace DistTool {

/**
 * Statistical utilities for SOAP descriptor analysis.
 *
 * Implements the metrics described in the paper:
 *   - Mean descriptor vector  q̄(T)  from a thermalized reference sample.
 *   - Covariance matrix  Σ  and its inverse.
 *   - Euclidean distance  d^i = ||q̃^i − q̄(T)||  (Eq. 5).
 *   - Mahalanobis distance  sqrt( (q̃^i−q̄)ᵀ Σ⁻¹ (q̃^i−q̄) ).
 *   - Gaussian defect probability  P(q̃^i | T)  (Eq. 6).
 */
class Statistics {
public:
    // Arithmetic mean of a set of vectors.
    static std::vector<double> mean(
        const std::vector<std::vector<double>>& data);

    // Unbiased sample covariance matrix (Bessel-corrected, n−1 denominator).
    static Eigen::MatrixXd covariance(
        const std::vector<std::vector<double>>& data,
        const std::vector<double>&              mean_vec);

    // Euclidean distance  ||a − b||.
    static double euclidean(const std::vector<double>& a,
                             const std::vector<double>& b);

    // Mahalanobis distance  sqrt( Δᵀ Σ⁻¹ Δ ),  Δ = x − μ.
    static double mahalanobis(const std::vector<double>& x,
                               const std::vector<double>& mu,
                               const Eigen::MatrixXd&     inv_cov);

    // Gaussian-model probability of being a lattice atom:
    //   P(d | T) = exp( −(d − mean_dist)² / (2 · var_dist) )
    // Centered on mean_dist so that a typical reference atom gives P ≈ 1.
    static double latticeProbability(double dist, double mean_dist, double var_dist);

    // Mean and variance of ||q̃^i − q̄|| over the reference sample.
    static void distanceStats(
        const std::vector<std::vector<double>>& dvs,
        const std::vector<double>&              mean_dv,
        double& mean_dist, double& var_dist);
};

} // namespace DistTool
