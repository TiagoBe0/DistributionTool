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

    // Mean and variance of ||q̃^i − q̄|| over the reference sample.
    static void distanceStats(
        const std::vector<std::vector<double>>& dvs,
        const std::vector<double>&              mean_dv,
        double& mean_dist, double& var_dist);
};

} // namespace DistTool
