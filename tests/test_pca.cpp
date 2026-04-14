#include "test_harness.h"
#include "AtomData.h"
#include "PCA.h"
#include <cmath>
#include <vector>

using namespace DistTool;

static std::vector<Atom> make_pca_atoms(const std::vector<std::vector<double>>& dvs) {
    std::vector<Atom> atoms(dvs.size());
    for (int i = 0; i < (int)dvs.size(); ++i) atoms[i].dv = dvs[i];
    return atoms;
}

// ── fit ───────────────────────────────────────────────────────────────────────

TEST_CASE("PCA::fit — throws on empty data") {
    PCA pca;
    CHECK_THROWS(pca.fit(std::vector<Atom>{}));
}

TEST_CASE("PCA::fit — throws if fewer than 2 samples") {
    PCA pca;
    Atom a; a.dv = {1.0, 2.0};
    CHECK_THROWS(pca.fit({a}));
}

TEST_CASE("PCA::fit — explained variance ratios sum to 1") {
    std::vector<std::vector<double>> dvs;
    for (int i = 0; i < 50; ++i)
        dvs.push_back({double(i), double(i) * 0.5, double(i) * 0.1});
    PCA pca;
    pca.fit(make_pca_atoms(dvs));
    const auto& evr = pca.explainedVarianceRatio();
    double sum = 0.0;
    for (double v : evr) sum += v;
    APPROX_EQ(sum, 1.0, 1e-10);
}

TEST_CASE("PCA::fit — 1D data: first component explains all variance") {
    // All DVs lie exactly on a line → PC1 should explain essentially 100%.
    std::vector<std::vector<double>> dvs;
    for (int i = 0; i < 30; ++i)
        dvs.push_back({double(i), 2.0 * double(i), 0.5 * double(i)});
    PCA pca;
    pca.fit(make_pca_atoms(dvs));
    const auto& evr = pca.explainedVarianceRatio();
    REQUIRE(evr[0] > 0.999);
}

// ── transform ─────────────────────────────────────────────────────────────────

TEST_CASE("PCA::transform — throws if not fitted") {
    PCA pca;
    Atom a; a.dv = {1.0, 2.0};
    CHECK_THROWS(pca.transform({a}, 1));
}

TEST_CASE("PCA::transform — projection of mean atom is zero") {
    // After centering, the mean DV projects to zero on every principal component.
    const int N = 30, D = 4;
    std::vector<std::vector<double>> dvs(N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
            dvs[i].push_back(double(i * (j + 1)));

    auto atoms = make_pca_atoms(dvs);
    PCA pca;
    pca.fit(atoms);

    // Compute mean DV manually
    std::vector<double> mean_dv(D, 0.0);
    for (const auto& a : atoms)
        for (int j = 0; j < D; ++j) mean_dv[j] += a.dv[j];
    for (auto& v : mean_dv) v /= N;

    Atom mean_atom; mean_atom.dv = mean_dv;
    auto proj = pca.transform({mean_atom}, 2);
    APPROX_EQ(proj[0][0], 0.0, 1e-9);
    APPROX_EQ(proj[0][1], 0.0, 1e-9);
}

TEST_CASE("PCA::transform — distance between same atom is zero") {
    std::vector<std::vector<double>> dvs;
    for (int i = 0; i < 20; ++i)
        dvs.push_back({double(i), double(i*i), 1.0 - double(i) * 0.05});
    auto atoms = make_pca_atoms(dvs);
    PCA pca;
    pca.fit(atoms);
    auto proj1 = pca.transform(atoms, 2);
    auto proj2 = pca.transform(atoms, 2);
    for (int i = 0; i < (int)atoms.size(); ++i) {
        APPROX_EQ(proj1[i][0], proj2[i][0], 1e-12);
        APPROX_EQ(proj1[i][1], proj2[i][1], 1e-12);
    }
}
