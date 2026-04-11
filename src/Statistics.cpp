#include "Statistics.h"
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
    for (const auto& v : data)
        for (int i = 0; i < d; ++i) m[i] += v[i];
    for (auto& x : m) x /= n;
    return m;
}

Eigen::MatrixXd Statistics::covariance(
    const std::vector<std::vector<double>>& data,
    const std::vector<double>&              mean_vec)
{
    if (data.size() < 2) throw std::invalid_argument("Need at least 2 samples for covariance");
    const int n = static_cast<int>(data.size());
    const int d = static_cast<int>(mean_vec.size());

    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(d, d);
    for (const auto& v : data) {
        Eigen::Map<const Eigen::VectorXd> xv(v.data(), d);
        Eigen::Map<const Eigen::VectorXd> mu(mean_vec.data(), d);
        Eigen::VectorXd diff = xv - mu;
        cov += diff * diff.transpose();
    }
    cov /= (n - 1);
    return cov;
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

double Statistics::mahalanobis(const std::vector<double>& x,
                                const std::vector<double>& mu,
                                const Eigen::MatrixXd&     inv_cov)
{
    const int d = static_cast<int>(x.size());
    Eigen::VectorXd diff(d);
    for (int i = 0; i < d; ++i) diff(i) = x[i] - mu[i];
    const double val = diff.transpose() * inv_cov * diff;
    return std::sqrt(std::max(0.0, val));
}

double Statistics::latticeProbability(double dist, double mean_dist, double var_dist) {
    // P(d|T) = exp( -(d - mean_dist)² / (2·var_dist) )
    // Centred on the reference mean so typical lattice atoms give P ≈ 1.
    if (var_dist < 1e-20) return (std::abs(dist - mean_dist) < 1e-10) ? 1.0 : 0.0;
    const double delta = dist - mean_dist;
    return std::exp(-delta * delta / (2.0 * var_dist));
}

void Statistics::distanceStats(
    const std::vector<std::vector<double>>& dvs,
    const std::vector<double>&              mean_dv,
    double& mean_dist, double& var_dist)
{
    const int n = static_cast<int>(dvs.size());
    if (n == 0) { mean_dist = var_dist = 0.0; return; }

    std::vector<double> dists(n);
    for (int i = 0; i < n; ++i)
        dists[i] = euclidean(dvs[i], mean_dv);

    mean_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / n;

    double var = 0.0;
    for (double d : dists) {
        const double delta = d - mean_dist;
        var += delta * delta;
    }
    var_dist = var / n;   // variance of distances over the reference sample
}

} // namespace DistTool
