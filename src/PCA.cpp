#include "PCA.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <stdexcept>

namespace DistTool {

void PCA::fit(const std::vector<Atom>& atoms) {
    if (atoms.empty()) throw std::invalid_argument("PCA::fit — empty data");

    const int n = static_cast<int>(atoms.size());
    const int d = static_cast<int>(atoms[0].dv.size());

    if (n < 2) throw std::invalid_argument("PCA::fit — need at least 2 samples");

    // ── Center the data ───────────────────────────────────────────────────────
    mean_ = Eigen::VectorXd::Zero(d);
    for (const auto& a : atoms)
        for (int j = 0; j < d; ++j) mean_(j) += a.dv[j];
    mean_ /= n;

    // ── Build centred data matrix X  [n × d] ─────────────────────────────────
    Eigen::MatrixXd X(n, d);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < d; ++j)
            X(i, j) = atoms[i].dv[j] - mean_(j);

    // ── Covariance matrix  C = XᵀX / (n−1),  [d × d] ────────────────────────
    Eigen::MatrixXd C = (X.transpose() * X) / static_cast<double>(n - 1);

    // ── Eigen-decomposition (symmetric → real eigenvalues) ───────────────────
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(C);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("PCA::fit — eigendecomposition failed");

    // Eigen returns eigenvalues in ascending order → reverse to descending
    const int nc = std::min(n - 1, d);   // at most n-1 non-zero components
    evals_ = solver.eigenvalues().tail(nc).reverse();
    // evecs_ rows = principal components (each row is one component)
    evecs_ = solver.eigenvectors().rightCols(nc).rowwise().reverse().transpose();

    // ── Explained variance ratio ──────────────────────────────────────────────
    const double total_var = solver.eigenvalues().sum();
    evr_.resize(nc);
    for (int k = 0; k < nc; ++k)
        evr_[k] = (total_var > 0.0) ? (evals_(k) / total_var) : 0.0;

    fitted_ = true;
}

std::vector<std::vector<double>> PCA::transform(
    const std::vector<Atom>& atoms,
    int n_components) const
{
    if (!fitted_) throw std::runtime_error("PCA::transform — not fitted yet");

    const int n  = static_cast<int>(atoms.size());
    const int d  = static_cast<int>(mean_.size());
    const int nc = std::min(n_components, static_cast<int>(evecs_.rows()));

    std::vector<std::vector<double>> result(n, std::vector<double>(nc, 0.0));

    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd x(d);
        for (int j = 0; j < d; ++j) x(j) = atoms[i].dv[j] - mean_(j);

        for (int k = 0; k < nc; ++k)
            result[i][k] = evecs_.row(k).dot(x);
    }

    return result;
}

} // namespace DistTool
