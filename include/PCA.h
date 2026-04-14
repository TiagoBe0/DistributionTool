#pragma once
#include "AtomData.h"
#include <Eigen/Dense>
#include <vector>

namespace DistTool {

/**
 * Principal Component Analysis (PCA) on descriptor vectors.
 *
 * Used (§3.2, Fig. 7) to visualise the distribution of DVs and to
 * discover unexpected defect types (type-a) not covered by the reference set.
 *
 * Steps:
 *   1. fit()        — center data, compute covariance, eigen-decompose.
 *   2. transform()  — project data onto the top-k principal components.
 *
 * Eigen-decomposition is performed with Eigen's SelfAdjointEigenSolver
 * (symmetric matrix → real eigenvalues, orthonormal eigenvectors).
 */
class PCA {
public:
    // Fit PCA to atoms — reads atom.dv directly, no intermediate copy.
    void fit(const std::vector<Atom>& atoms);

    // Project atoms onto the first `n_components` principal axes.
    // Returns matrix of shape [n_atoms][n_components].
    std::vector<std::vector<double>> transform(
        const std::vector<Atom>& atoms,
        int n_components = 2) const;

    // Explained variance fraction for each component (descending).
    const std::vector<double>& explainedVarianceRatio() const { return evr_; }

    // Eigenvalues in descending order.
    const Eigen::VectorXd& eigenvalues() const { return evals_; }

    // Principal components: row k = k-th component, shape [n_comp][n_features].
    const Eigen::MatrixXd& components() const { return evecs_; }

    bool fitted() const { return fitted_; }

private:
    Eigen::VectorXd    mean_;
    Eigen::MatrixXd    evecs_;    // rows = components (descending variance)
    Eigen::VectorXd    evals_;    // eigenvalues, descending
    std::vector<double> evr_;     // explained variance ratios
    bool               fitted_ = false;
};

} // namespace DistTool
