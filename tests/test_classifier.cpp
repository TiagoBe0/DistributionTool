#include "test_harness.h"
#include "AtomData.h"
#include "DefectClassifier.h"
#include <array>
#include <cmath>
#include <limits>
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

    auto frame = make_frame_with_dvs(std::vector<std::vector<double>>(200, dv));

    DefectClassifier clf(0.15);
    clf.buildReference(frame.atoms);
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

    auto ref_frame = make_frame_with_dvs(std::vector<std::vector<double>>(200, dv_ref));

    std::vector<std::vector<double>> all_dvs(200, dv_ref);
    all_dvs.push_back(dv_def);
    auto dmg_frame = make_frame_with_dvs(all_dvs);

    DefectClassifier clf(0.15);
    clf.buildReference(ref_frame.atoms);
    clf.classify(dmg_frame);

    // Last atom is the outlier
    const auto& outlier = dmg_frame.atoms.back();
    REQUIRE(outlier.defect_type != DefectType::Lattice);
    REQUIRE(outlier.dist_to_ref > 0.15);
    // defect_prob should be close to 1 for a very distorted atom
    REQUIRE(outlier.defect_prob > 0.9);
}

// Build a frame from explicit (type, dv) pairs.
static Frame make_typed_frame(const std::vector<int>& types,
                              const std::vector<std::vector<double>>& dvs) {
    Frame f;
    f.box = make_box(50.0);
    f.atoms.resize(dvs.size());
    for (int i = 0; i < (int)dvs.size(); ++i) {
        f.atoms[i].id   = i + 1;
        f.atoms[i].type = types[i];
        f.atoms[i].x    = i * 2.0;
        f.atoms[i].dv   = dvs[i];
    }
    return f;
}

// ── per-species reference ─────────────────────────────────────────────────────

TEST_CASE("DefectClassifier — per-species builds one reference per type") {
    // Type 1 atoms ≈ (1,0); type 2 atoms ≈ (0,1). A single global reference
    // would sit at (0.5,0.5), so every atom is ~0.707 from it → all flagged.
    // Per-species references sit exactly on each cluster → distance 0 → Lattice.
    std::vector<int> types;
    std::vector<std::vector<double>> dvs;
    for (int i = 0; i < 100; ++i) { types.push_back(1); dvs.push_back({1.0, 0.0}); }
    for (int i = 0; i < 100; ++i) { types.push_back(2); dvs.push_back({0.0, 1.0}); }
    auto frame = make_typed_frame(types, dvs);

    // Global reference: everything looks like a defect.
    {
        auto g = frame;
        DefectClassifier clf(0.15);
        clf.buildReference(g.atoms);
        clf.classify(g);
        int non_lattice = 0;
        for (const auto& a : g.atoms)
            if (a.defect_type != DefectType::Lattice) ++non_lattice;
        REQUIRE(non_lattice == 200);
    }

    // Per-species reference: every atom matches its species mean → Lattice.
    {
        auto p = frame;
        DefectClassifier clf(0.15);
        clf.setPerSpecies(true);
        clf.buildReference(p.atoms);
        REQUIRE(clf.referencesByType().size() == 2);
        clf.classify(p);
        for (const auto& a : p.atoms) {
            REQUIRE(a.defect_type == DefectType::Lattice);
            APPROX_EQ(a.dist_to_ref, 0.0, 1e-12);
        }
    }
}

TEST_CASE("DefectClassifier — under-populated species folds into global reference") {
    // 200 atoms of type 1 plus a single type-4 atom (like a tagged PKA). The
    // lone type-4 atom must NOT get its own (degenerate) per-species reference.
    std::vector<int> types(200, 1);
    std::vector<std::vector<double>> dvs(200, {1.0, 0.0});
    types.push_back(4);
    dvs.push_back({1.0, 0.0});
    auto ref = make_typed_frame(types, dvs);

    DefectClassifier clf(0.15);
    clf.setPerSpecies(true);
    clf.setThresholdPercentile(0.99);
    clf.buildReference(ref.atoms);

    // Only type 1 is populous enough for a per-species reference.
    REQUIRE(clf.referencesByType().size() == 1);
    REQUIRE(clf.referencesByType().count(1) == 1);
    REQUIRE(clf.referencesByType().count(4) == 0);
}

TEST_CASE("DefectClassifier — unseen type falls back to global reference") {
    // Reference has only type 1; classify an atom of an unseen type 9.
    // It must use the global reference (no crash, deterministic distance).
    std::vector<int> types(50, 1);
    std::vector<std::vector<double>> dvs(50, {1.0, 0.0});
    auto ref = make_typed_frame(types, dvs);

    DefectClassifier clf(0.15);
    clf.setPerSpecies(true);
    clf.buildReference(ref.atoms);

    Frame dmg = ref;
    dmg.atoms.push_back({});
    dmg.atoms.back().id = 999;
    dmg.atoms.back().type = 9;        // unseen species
    dmg.atoms.back().dv = {1.0, 0.0}; // identical to global mean → dist 0
    clf.classify(dmg);
    APPROX_EQ(dmg.atoms.back().dist_to_ref, 0.0, 1e-12);
    REQUIRE(dmg.atoms.back().defect_type == DefectType::Lattice);
}

// ── empirical defect probability ──────────────────────────────────────────────

TEST_CASE("DefectClassifier — empirical defect_prob is a percentile rank") {
    // Reference DVs along axis 0 at {0,1,2,3,4}; mean = (2,0).
    // Reference distances = {2,1,0,1,2} → sorted {0,1,1,2,2}.
    auto ref = make_typed_frame({1,1,1,1,1},
        {{0.0,0.0},{1.0,0.0},{2.0,0.0},{3.0,0.0},{4.0,0.0}});

    DefectClassifier clf(0.15);
    clf.setEmpiricalProb(true);
    clf.buildReference(ref.atoms);

    // Damaged atoms with distances 0, 1, 2 from the mean (2,0).
    auto dmg = make_typed_frame({1,1,1},
        {{2.0,0.0},{1.0,0.0},{0.0,0.0}});
    clf.classify(dmg);

    APPROX_EQ(dmg.atoms[0].dist_to_ref, 0.0, 1e-12);
    APPROX_EQ(dmg.atoms[1].dist_to_ref, 1.0, 1e-12);
    APPROX_EQ(dmg.atoms[2].dist_to_ref, 2.0, 1e-12);
    // empiricalCdf over {0,1,1,2,2}: F(0)=1/5, F(1)=3/5, F(2)=5/5.
    APPROX_EQ(dmg.atoms[0].defect_prob, 0.2, 1e-12);
    APPROX_EQ(dmg.atoms[1].defect_prob, 0.6, 1e-12);
    APPROX_EQ(dmg.atoms[2].defect_prob, 1.0, 1e-12);
}

// ── percentile-derived threshold ──────────────────────────────────────────────

TEST_CASE("DefectClassifier — threshold-percentile sets the global threshold") {
    // Reference distances sorted {0,1,1,2,2}; median (p=0.5) = 1.0.
    auto ref = make_typed_frame({1,1,1,1,1},
        {{0.0,0.0},{1.0,0.0},{2.0,0.0},{3.0,0.0},{4.0,0.0}});

    DefectClassifier clf(0.15);              // fixed threshold ignored once pct set
    clf.setThresholdPercentile(0.5);
    clf.buildReference(ref.atoms);
    APPROX_EQ(clf.reference().threshold, 1.0, 1e-12);

    // d < 1 → Lattice;  d >= 1 → non-Lattice.
    auto dmg = make_typed_frame({1,1}, {{2.0,0.0},{0.0,0.0}});  // d=0, d=2
    clf.classify(dmg);
    REQUIRE(dmg.atoms[0].defect_type == DefectType::Lattice);     // d=0 < 1
    REQUIRE(dmg.atoms[1].defect_type != DefectType::Lattice);     // d=2 >= 1
}

TEST_CASE("DefectClassifier — per-species percentile thresholds are independent") {
    // Type 1 distances spread to ~1; type 2 distances spread to ~10.
    // Repeat each 5-value pattern 20× → 100 atoms/species (above the minimum)
    // while preserving the median distance (1.0 for type 1, 10.0 for type 2).
    std::vector<int> types;
    std::vector<std::vector<double>> dvs;
    for (int rep = 0; rep < 20; ++rep) {
        for (double v : {0.0,1.0,2.0,3.0,4.0})     { types.push_back(1); dvs.push_back({v,0.0}); }
        for (double v : {0.0,10.0,20.0,30.0,40.0}) { types.push_back(2); dvs.push_back({v,0.0}); }
    }
    auto ref = make_typed_frame(types, dvs);

    DefectClassifier clf(0.15);
    clf.setPerSpecies(true);
    clf.setThresholdPercentile(0.5);
    clf.buildReference(ref.atoms);

    const auto& by_type = clf.referencesByType();
    REQUIRE(by_type.size() == 2);
    APPROX_EQ(by_type.at(1).threshold,  1.0, 1e-12);   // median dist for type 1
    APPROX_EQ(by_type.at(2).threshold, 10.0, 1e-12);   // median dist for type 2
}

TEST_CASE("DefectClassifier::buildReference — throws on empty DVs") {
    DefectClassifier clf(0.15);
    Frame f;
    f.box = make_box(10.0);
    Atom a; a.id = 1; a.type = 1; a.x = 1.0; a.y = 1.0; a.z = 1.0;
    // dv is default-empty
    f.atoms.push_back(a);
    CHECK_THROWS(clf.buildReference(f.atoms));
}

// ── findVacanciesGrid ─────────────────────────────────────────────────────────

TEST_CASE("findVacanciesGrid — throws on non-positive spacing or threshold") {
    DefectClassifier clf(0.15);
    Frame f; f.box = make_box(10.0);
    CHECK_THROWS(clf.findVacanciesGrid(f, 0.0, 1.0));
    CHECK_THROWS(clf.findVacanciesGrid(f, 1.0, 0.0));
    CHECK_THROWS(clf.findVacanciesGrid(f, -1.0, 1.0));
}

TEST_CASE("findVacanciesGrid — empty frame: no spurious overflow vacancies") {
    // No atoms → CellList returns DBL_MAX, which is filtered out so the result
    // is empty (instead of emitting d_near values like sqrt(DBL_MAX) ≈ 1.3e154
    // that pollute downstream CSVs and clustering).
    DefectClassifier clf(0.15);
    Frame f; f.box = make_box(10.0);
    auto vac = clf.findVacanciesGrid(f, 2.5, 1.0);
    REQUIRE(vac.empty());
}

TEST_CASE("findVacanciesGrid — dense atom packing: no vacancies") {
    // Place atoms on a 1 Å grid in a 5×5×5 box; threshold = 0.9 Å.
    // Every grid point (spacing 1 Å) is within 0.9 Å of at least one atom.
    DefectClassifier clf(0.15);
    Frame f; f.box = make_box(5.0);
    int id = 1;
    for (int ix = 0; ix < 5; ++ix)
    for (int iy = 0; iy < 5; ++iy)
    for (int iz = 0; iz < 5; ++iz) {
        Atom a; a.id = id++; a.type = 1;
        a.x = 0.5 + ix; a.y = 0.5 + iy; a.z = 0.5 + iz;
        f.atoms.push_back(a);
    }
    auto vac = clf.findVacanciesGrid(f, 1.0, 0.6);
    REQUIRE(vac.empty());
}

TEST_CASE("findVacanciesGrid — single atom: vacancy points are far from atom") {
    DefectClassifier clf(0.15);
    Frame f; f.box = make_box(10.0);
    Atom a; a.id=1; a.type=1; a.x=5.0; a.y=5.0; a.z=5.0;
    f.atoms.push_back(a);

    const double threshold = 3.0;
    auto vac = clf.findVacanciesGrid(f, 1.0, threshold);

    for (const auto& v : vac) {
        double dx = v.pos[0]-5.0, dy = v.pos[1]-5.0, dz = v.pos[2]-5.0;
        REQUIRE(std::sqrt(dx*dx + dy*dy + dz*dz) > threshold);
        REQUIRE(v.d_near > threshold);
    }
}
