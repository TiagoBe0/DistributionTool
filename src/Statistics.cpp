#include "Statistics.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace DistTool {

std::vector<double> Statistics::mean(
    const std::vector<std::vector<double>>& data)
{
    if (data.empty()) return {};
    const int d = static_cast<int>(data[0].size());
    const int n = static_cast<int>(data.size());

    std::vector<double> m(d, 0.0);
    for (const auto& v : data) {
        if ((int)v.size() != d)
            throw std::invalid_argument(
                "Statistics::mean: inconsistent vector sizes (" +
                std::to_string(v.size()) + " vs " + std::to_string(d) + ")");
        for (int i = 0; i < d; ++i) m[i] += v[i];
    }
    for (auto& x : m) x /= n;
    return m;
}

double Statistics::euclidean(const std::vector<double>& a,
                              const std::vector<double>& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vector size mismatch in euclidean()");
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double d = a[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

double Statistics::chiProbability(double d, double k, double sigma) {
    // Chi-distribution probability, normalised to 1 at the mode.
    //
    // PDF shape:  f(d) ∝ d^(k-2) · exp(-d²/(2σ²))
    // Mode at d_peak = σ·√(k-2)  for k > 2
    //        d_peak = 0           for k ≤ 2  (Rayleigh / exponential shape)
    //
    // Normalised form: P(d) = f(d) / f(d_peak) so that P ∈ (0, 1]
    // and P = 1 exactly at the reference mean environment.
    if (sigma < 1e-20) return (d < 1e-10) ? 1.0 : 0.0;
    if (d < 0.0) return 0.0;

    if (k <= 2.0) {
        // Mode at 0 → monotone decreasing: P = exp(-d²/(2σ²))
        return std::exp(-d * d / (2.0 * sigma * sigma));
    }

    // Mode at d_peak = σ·√(k-2)
    const double d_peak = sigma * std::sqrt(k - 2.0);
    // log P(d) - log P(d_peak) = (k-2)·log(d/d_peak) - (d²-d_peak²)/(2σ²)
    if (d < 1e-20) return 0.0;   // d^(k-2) → 0 for k > 2 when d→0
    const double log_ratio = (k - 2.0) * std::log(d / d_peak)
                           - (d * d - d_peak * d_peak) / (2.0 * sigma * sigma);
    return std::exp(log_ratio);
}

void Statistics::fitChiParams(const std::vector<double>& dists,
                               double& k_out, double& sigma_out)
{
    // Method of moments for the chi distribution.
    // For d ~ chi(k, σ):  E[d²] = k·σ²,  Var(d²) = 2k·σ⁴
    // → k = 2·(E[d²])² / Var(d²),   σ² = E[d²] / k
    const int n = static_cast<int>(dists.size());
    if (n < 2) { k_out = 2.0; sigma_out = 1.0; return; }

    double mean_d2 = 0.0, mean_d4 = 0.0;
    for (double d : dists) {
        const double d2 = d * d;
        mean_d2 += d2;
        mean_d4 += d2 * d2;
    }
    mean_d2 /= n;
    mean_d4 /= n;

    const double var_d2 = mean_d4 - mean_d2 * mean_d2;
    if (var_d2 < 1e-20 || mean_d2 < 1e-20) {
        // Degenerate: fall back to Rayleigh (k=2)
        k_out     = 2.0;
        sigma_out = std::sqrt(mean_d2 / 2.0);
        return;
    }
    // Cap k to avoid catastrophic cancellation in chiProbability for very
    // sharply-peaked distributions (e.g. near-perfect crystals at low T).
    // In practice k ≈ effective DV dimensionality; 2000 is generous.
    constexpr double k_max = 2000.0;
    k_out     = std::min(k_max, std::max(2.0, 2.0 * mean_d2 * mean_d2 / var_d2));
    sigma_out = std::sqrt(mean_d2 / k_out);
}

double Statistics::empiricalCdf(const std::vector<double>& sorted_asc, double d)
{
    const size_t n = sorted_asc.size();
    if (n == 0) return 0.0;
    // Count of reference samples ≤ d, divided by n.
    const auto it = std::upper_bound(sorted_asc.begin(), sorted_asc.end(), d);
    return static_cast<double>(it - sorted_asc.begin()) / static_cast<double>(n);
}

double Statistics::percentile(const std::vector<double>& sorted_asc, double p)
{
    const size_t n = sorted_asc.size();
    if (n == 0) return 0.0;
    if (n == 1) return sorted_asc[0];
    if (p <= 0.0) return sorted_asc.front();
    if (p >= 1.0) return sorted_asc.back();
    // Linear interpolation between order statistics (NumPy 'linear' method).
    const double pos = p * (n - 1);
    const size_t lo  = static_cast<size_t>(pos);
    const double frac = pos - lo;
    if (lo + 1 >= n) return sorted_asc.back();
    return sorted_asc[lo] + frac * (sorted_asc[lo + 1] - sorted_asc[lo]);
}

} // namespace DistTool
