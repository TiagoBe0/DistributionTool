#include "test_harness.h"
#include "AtomData.h"
#include "DefectClassifier.h"
#include <array>
#include <cmath>
#include <vector>

using namespace DistTool;

// ── clusterVacancyPoints ──────────────────────────────────────────────────────

// Helper: non-periodic box [0, L]^3
static SimBox make_box(double L) {
    SimBox b;
    b.xb = {0, L};  b.yb = {0, L};  b.zb = {0, L};
    b.periodic[0] = b.periodic[1] = b.periodic[2] = false;
    return b;
}

TEST_CASE("clusterVacancyPoints — empty input gives empty output") {
    DefectClassifier clf(0.15);
    auto clusters = clf.clusterVacancyPoints({}, 1.0, make_box(20.0));
    REQUIRE(clusters.empty());
}

TEST_CASE("clusterVacancyPoints — single point gives one cluster") {
    DefectClassifier clf(0.15);
    std::vector<VacancyPoint> pts = {{{5.0, 5.0, 5.0}, 2.5}};
    auto clusters = clf.clusterVacancyPoints(pts, 2.0, make_box(20.0));
    REQUIRE(clusters.size() == 1);
    REQUIRE(clusters[0].n_pts == 1);
    APPROX_EQ(clusters[0].d_near_max, 2.5, 1e-12);
}

TEST_CASE("clusterVacancyPoints — two close points merge into one cluster") {
    DefectClassifier clf(0.15);
    std::vector<VacancyPoint> pts = {
        {{5.0, 5.0, 5.0}, 2.5},
        {{5.5, 5.0, 5.0}, 2.3}   // distance = 0.5 < r_cluster=1.0
    };
    auto clusters = clf.clusterVacancyPoints(pts, 1.0, make_box(20.0));
    REQUIRE(clusters.size() == 1);
    REQUIRE(clusters[0].n_pts == 2);
    // Cluster seeded from deepest point (d_near=2.5)
    APPROX_EQ(clusters[0].d_near_max, 2.5, 1e-12);
}

TEST_CASE("clusterVacancyPoints — two distant points give two clusters") {
    DefectClassifier clf(0.15);
    std::vector<VacancyPoint> pts = {
        {{ 2.0,  2.0,  2.0}, 2.5},
        {{18.0, 18.0, 18.0}, 2.3}   // distance >> r_cluster=1.0
    };
    auto clusters = clf.clusterVacancyPoints(pts, 1.0, make_box(20.0));
    REQUIRE(clusters.size() == 2);
    REQUIRE(clusters[0].n_pts == 1);
    REQUIRE(clusters[1].n_pts == 1);
}

TEST_CASE("clusterVacancyPoints — cluster seeded from deepest point") {
    DefectClassifier clf(0.15);
    // Three points: A (deep) and B (shallow) are far apart; B and C are close.
    // Correct order: seed A first (d_near=3.0), then seed B/C cluster.
    std::vector<VacancyPoint> pts = {
        {{ 2.0, 2.0, 2.0}, 3.0},   // A — deepest
        {{15.0, 2.0, 2.0}, 1.5},   // B — shallow, close to C
        {{15.8, 2.0, 2.0}, 1.2}    // C — shallowest, close to B
    };
    auto clusters = clf.clusterVacancyPoints(pts, 2.0, make_box(20.0));
    REQUIRE(clusters.size() == 2);
    // First cluster must be the isolated deep point A
    APPROX_EQ(clusters[0].d_near_max, 3.0, 1e-12);
    REQUIRE(clusters[0].n_pts == 1);
    // Second cluster is B+C merged
    APPROX_EQ(clusters[1].d_near_max, 1.5, 1e-12);
    REQUIRE(clusters[1].n_pts == 2);
}

// ── DefectClassifier pipeline ─────────────────────────────────────────────────

// Build a minimal Frame where each atom already has its dv set.
static Frame make_frame_with_dvs(const std::vector<std::vector<double>>& dvs) {
    Frame f;
    f.box = make_box(50.0);
    f.atoms.resize(dvs.size());
    for (int i = 0; i < (int)dvs.size(); ++i) {
        f.atoms[i].id   = i + 1;
        f.atoms[i].type = 1;
        f.atoms[i].x    = i * 2.0;
        f.atoms[i].y    = 0.0;
        f.atoms[i].z    = 0.0;
        f.atoms[i].dv   = dvs[i];
    }
    return f;
}

TEST_CASE("DefectClassifier::classify — identical DVs all become Lattice") {
    // If every atom has the same DV, dist_to_ref = 0 for all → all Lattice.
    const int D = 10;
    std::vector<double> dv(D, 0.0);
    dv[0] = 1.0;  // unit vector along first axis

    std::vector<std::vector<double>> ref_dvs(200, dv);
    auto frame = make_frame_with_dvs(ref_dvs);

    DefectClassifier clf(0.15);
    clf.buildReference(ref_dvs);
    clf.classify(frame);

    for (const auto& a : frame.atoms) {
        REQUIRE(a.defect_type == DefectType::Lattice);
        APPROX_EQ(a.dist_to_ref, 0.0, 1e-12);
    }
}

TEST_CASE("DefectClassifier::classify — very different DV flagged as non-Lattice") {
    const int D = 10;
    // Reference: unit vector along axis 0
    std::vector<double> dv_ref(D, 0.0);  dv_ref[0] = 1.0;
    // Defect: unit vector along axis 1 — distance = sqrt(2) >> threshold
    std::vector<double> dv_def(D, 0.0);  dv_def[1] = 1.0;

    std::vector<std::vector<double>> ref_dvs(200, dv_ref);
    std::vector<std::vector<double>> all_dvs = ref_dvs;
    all_dvs.push_back(dv_def);

    auto frame = make_frame_with_dvs(all_dvs);

    DefectClassifier clf(0.15);
    clf.buildReference(ref_dvs);
    clf.classify(frame);

    // Last atom is the outlier
    const auto& outlier = frame.atoms.back();
    REQUIRE(outlier.defect_type != DefectType::Lattice);
    REQUIRE(outlier.dist_to_ref > 0.15);
    // defect_prob should be close to 1 for a very distorted atom
    REQUIRE(outlier.defect_prob > 0.9);
}
