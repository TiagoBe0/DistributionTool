#include "test_harness.h"
#include "AtomData.h"
#include "HybridVacancyDetector.h"
#include "WignerSeitz.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace DistTool;

// ── Helpers ───────────────────────────────────────────────────────────────────

// Build a BCC pristine lattice of nc×nc×nc unit cells with lattice constant a.
// Returns a Frame with periodic boundary conditions.
static Frame make_bcc(int nc, double a) {
    Frame f;
    f.box.xb = {0.0, nc * a};
    f.box.yb = {0.0, nc * a};
    f.box.zb = {0.0, nc * a};
    f.box.periodic[0] = f.box.periodic[1] = f.box.periodic[2] = true;

    int id = 1;
    for (int ix = 0; ix < nc; ++ix)
    for (int iy = 0; iy < nc; ++iy)
    for (int iz = 0; iz < nc; ++iz) {
        Atom corner;
        corner.id = id++;
        corner.type = 1;
        corner.x = ix * a;
        corner.y = iy * a;
        corner.z = iz * a;
        f.atoms.push_back(corner);

        Atom center;
        center.id = id++;
        center.type = 1;
        center.x = (ix + 0.5) * a;
        center.y = (iy + 0.5) * a;
        center.z = (iz + 0.5) * a;
        f.atoms.push_back(center);
    }
    return f;
}

// Add Gaussian thermal noise to every atom (σ in Å).
static void thermalize(Frame& f, double sigma, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, sigma);
    for (auto& a : f.atoms) {
        a.x += nd(rng);
        a.y += nd(rng);
        a.z += nd(rng);
    }
}

// Remove atoms with the given indices (sorted-descending to keep indices valid).
static void remove_atoms(Frame& f, std::vector<int> indices) {
    std::sort(indices.begin(), indices.end(), std::greater<int>());
    for (int i : indices)
        f.atoms.erase(f.atoms.begin() + i);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("HybridParams::applyPreset — invalid name throws") {
    HybridParams p;
    CHECK_THROWS(p.applyPreset("nonexistent"));
}

TEST_CASE("HybridParams::applyPreset — 'ws' yields WS-only weights") {
    HybridParams p;
    p.applyPreset("ws");
    REQUIRE(p.w_ws == 1.0);
    REQUIRE(p.w_soft_ws == 0.0);
    REQUIRE(p.w_voronoi == 0.0);
    REQUIRE(p.w_density == 0.0);
    REQUIRE(p.w_soap == 0.0);
    REQUIRE(p.w_topology == 0.0);
    REQUIRE(p.w_transit == 0.0);
    REQUIRE(p.w_frenkel == 0.0);
}

TEST_CASE("Hybrid 'ws' preset — pristine lattice: 0 vacancies (matches WS)") {
    const int nc = 6;
    const double a = 2.87;
    Frame ref = make_bcc(nc, a);
    Frame dmg = ref;

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    REQUIRE(ws.vacancyCount() == 0);

    HybridParams p;
    p.applyPreset("ws");
    HybridVacancyDetector hvd(p);
    hvd.buildReference(ref, ws);
    auto vacs = hvd.detect(dmg, ws_res);
    REQUIRE(hvd.vacancyCount() == 0);
}

TEST_CASE("Hybrid 'ws' preset — 5 vacancies, T=0: count matches WS exact") {
    const int nc = 6;
    const double a = 2.87;
    Frame ref = make_bcc(nc, a);

    Frame dmg = ref;
    // Choose 5 atoms to delete (non-adjacent to avoid void merging).
    std::vector<int> to_remove = {10, 100, 200, 300, 400};
    remove_atoms(dmg, to_remove);

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    REQUIRE(ws.vacancyCount() == 5);

    HybridParams p;
    p.applyPreset("ws");
    HybridVacancyDetector hvd(p);
    hvd.buildReference(ref, ws);
    auto vacs = hvd.detect(dmg, ws_res);

    // Invariant: when only w_ws is non-zero, hybrid count must equal WS count.
    REQUIRE(hvd.vacancyCount() == ws.vacancyCount());
    REQUIRE(hvd.vacancyCount() == 5);
}

TEST_CASE("Hybrid 'robust' preset — 10 vacancies, T=0: recall == 1") {
    const int nc = 6;
    const double a = 2.87;
    Frame ref = make_bcc(nc, a);

    Frame dmg = ref;
    std::vector<int> to_remove = {5, 35, 70, 110, 150, 200, 260, 320, 380, 420};
    remove_atoms(dmg, to_remove);

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    REQUIRE(ws.vacancyCount() == 10);

    HybridParams p;
    p.applyPreset("robust");
    HybridVacancyDetector hvd(p);
    hvd.buildReference(ref, ws);
    auto vacs = hvd.detect(dmg, ws_res);

    // The robust consensus must include all true vacancies (no missed sites).
    // With no thermal noise and no transit atoms, all 10 should pass.
    REQUIRE(hvd.vacancyCount() >= 10);
    REQUIRE(hvd.wsAgreeCount() == 10);
}

TEST_CASE("Hybrid 'robust' preset — ballistic transit: WS overcounts, hybrid filters") {
    // Setup: pristine BCC, displace ONE atom by ~1.2*NN into a vacant interstitial-like
    // spot.  WS will count this as one vacancy (atom left its site) + one interstitial
    // (atom crowded another site).  The robust hybrid should filter the vacancy
    // via the transit penalty (atom is still near its original site).
    const int nc = 6;
    const double a = 2.87;
    const double nn = a * std::sqrt(3.0) / 2.0;     // BCC 1NN ≈ 2.484 Å

    Frame ref = make_bcc(nc, a);
    Frame dmg = ref;

    // Displace atom 50 by 1.15·nn along +x: just past 0.5·nn so WS sees it as
    // belonging to a different cell, but still inside the transit band.
    const int target = 50;
    dmg.atoms[target].x += 1.15 * nn;

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);

    // WS should now report at least one vacancy.
    REQUIRE(ws.vacancyCount() >= 1);

    // ws preset: hybrid must agree with WS exactly.
    {
        HybridParams p;
        p.applyPreset("ws");
        HybridVacancyDetector hvd(p);
        hvd.buildReference(ref, ws);
        hvd.detect(dmg, ws_res);
        REQUIRE(hvd.vacancyCount() == ws.vacancyCount());
    }

    // robust preset: the displaced atom is at ~1.15·nn from the vacated site →
    // inside the transit band [0.5·nn, 1.5·nn], so transit_penalty should fire
    // and the candidate should be filtered out.
    {
        HybridParams p;
        p.applyPreset("robust");
        HybridVacancyDetector hvd(p);
        hvd.buildReference(ref, ws);
        hvd.detect(dmg, ws_res);
        // The robust count should be strictly LESS than WS for this case.
        REQUIRE(hvd.vacancyCount() < ws.vacancyCount());
        REQUIRE(hvd.wsOnlyCount() >= 1);     // WS-only ≥ 1 (filtered the transit)
    }
}

TEST_CASE("HybridParams::applyPreset — 'relaxed' disables transit/Frenkel penalties") {
    HybridParams p;
    p.applyPreset("relaxed");
    REQUIRE(p.w_transit == 0.0);
    REQUIRE(p.w_frenkel == 0.0);
    // Positive multi-signal weights stay enabled (robust to thermal noise).
    REQUIRE(p.w_ws == 1.0);
    REQUIRE(p.w_soft_ws > 0.0);
}

TEST_CASE("Hybrid 'relaxed' preset — surviving Frenkel pair is kept, not filtered") {
    // Same geometry as the ballistic-transit case, but here we model a COOLED
    // final state: the displaced atom has settled at an interstitial site far
    // enough that it is a stable surviving defect. 'robust' filters the vacancy
    // (transit/Frenkel penalty); 'relaxed' must KEEP it and agree with WS.
    const int nc = 6;
    const double a = 2.87;
    const double nn = a * std::sqrt(3.0) / 2.0;

    Frame ref = make_bcc(nc, a);
    Frame dmg = ref;

    const int target = 50;
    dmg.atoms[target].x += 1.15 * nn;   // displaced into the transit band

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    REQUIRE(ws.vacancyCount() >= 1);

    // robust: filters the displaced-atom vacancy → count < WS.
    {
        HybridParams p;
        p.applyPreset("robust");
        HybridVacancyDetector hvd(p);
        hvd.buildReference(ref, ws);
        hvd.detect(dmg, ws_res);
        REQUIRE(hvd.vacancyCount() < ws.vacancyCount());
    }

    // relaxed: no penalties → every WS-vacant site is accepted, count == WS.
    {
        HybridParams p;
        p.applyPreset("relaxed");
        HybridVacancyDetector hvd(p);
        hvd.buildReference(ref, ws);
        hvd.detect(dmg, ws_res);
        REQUIRE(hvd.vacancyCount() == ws.vacancyCount());
        REQUIRE(hvd.wsOnlyCount() == 0);   // nothing filtered out
    }
}

TEST_CASE("HybridParams::applyPreset — 'survival' enables cloud signals, disables auto-accept") {
    HybridParams p;
    p.applyPreset("survival");
    REQUIRE(p.w_centrality > 0.0);
    REQUIRE(p.w_cloud_density > 0.0);
    REQUIRE(p.w_transit == 0.0);
    REQUIRE(p.w_frenkel == 0.0);
    REQUIRE(p.ws_auto_accept == false);

    // The historical presets must keep the cloud signals OFF and auto-accept ON
    // (preserves their WS-equivalence invariants).
    for (const char* name : {"ws", "robust", "relaxed", "sensitive"}) {
        HybridParams q;
        q.applyPreset(name);
        REQUIRE(q.w_centrality == 0.0);
        REQUIRE(q.w_cloud_density == 0.0);
        REQUIRE(q.ws_auto_accept == true);
    }
}

TEST_CASE("Hybrid cloud signals — clustered vacancies rank above an isolated one") {
    // 4 vacancies clustered near the box center + 1 isolated far away (> cloud
    // radius from the cluster). Core-shell expectation: the clustered sites get
    // higher cloud_centrality and cloud_density than the isolated one, and the
    // isolated one has the lowest consensus score under the 'survival' preset.
    const int nc = 8;                      // box 22.96 Å (half-box > cloud_radius)
    const double a = 2.87;
    Frame ref = make_bcc(nc, a);
    Frame dmg = ref;

    // Find 4 mutually-close atoms near the center and 1 far from all of them.
    const double cx = nc * a / 2.0;
    std::vector<int> cluster;
    for (int i = 0; i < (int)ref.atoms.size() && (int)cluster.size() < 4; ++i) {
        const auto& at = ref.atoms[i];
        const double r = std::sqrt((at.x-cx)*(at.x-cx) + (at.y-cx)*(at.y-cx)
                                   + (at.z-cx)*(at.z-cx));
        if (r < 1.1 * a) cluster.push_back(i);
    }
    REQUIRE((int)cluster.size() == 4);
    int isolated = -1;
    for (int i = 0; i < (int)ref.atoms.size(); ++i) {
        const auto& at = ref.atoms[i];
        // corner region: ≥ 11 Å from the center cluster in minimum image
        double dx = at.x - cx, dy = at.y - cx, dz = at.z - cx;
        const double L = nc * a;
        dx -= L * std::round(dx / L); dy -= L * std::round(dy / L);
        dz -= L * std::round(dz / L);
        if (std::sqrt(dx*dx + dy*dy + dz*dz) > 11.0) { isolated = i; break; }
    }
    REQUIRE(isolated >= 0);

    const double iso_x = ref.atoms[isolated].x;
    const double iso_y = ref.atoms[isolated].y;
    const double iso_z = ref.atoms[isolated].z;

    auto to_remove = cluster;
    to_remove.push_back(isolated);
    remove_atoms(dmg, to_remove);

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    REQUIRE(ws.vacancyCount() == 5);

    HybridParams p;
    p.applyPreset("survival");
    HybridVacancyDetector hvd(p);
    hvd.buildReference(ref, ws);
    auto vacs = hvd.detect(dmg, ws_res);

    // Locate the isolated candidate by position; compare against the others.
    const HybridVacancy* iso = nullptr;
    double min_clu_centrality = 2.0, min_clu_density = 2.0;
    for (const auto& v : vacs) {
        if (v.ref_site_idx < 0) continue;
        const double d = std::sqrt((v.pos[0]-iso_x)*(v.pos[0]-iso_x)
                                 + (v.pos[1]-iso_y)*(v.pos[1]-iso_y)
                                 + (v.pos[2]-iso_z)*(v.pos[2]-iso_z));
        if (d < 0.1) {
            iso = &v;
        } else {
            min_clu_centrality = std::min(min_clu_centrality, v.breakdown.cloud_centrality);
            min_clu_density    = std::min(min_clu_density,    v.breakdown.cloud_density);
        }
    }
    REQUIRE(iso != nullptr);
    REQUIRE(iso->breakdown.cloud_centrality < min_clu_centrality);
    REQUIRE(iso->breakdown.cloud_density    < min_clu_density);

    // The isolated vacancy must have the lowest consensus score of all 5.
    double min_other_score = 2.0;
    for (const auto& v : vacs) {
        if (v.ref_site_idx < 0 || &v == iso) continue;
        min_other_score = std::min(min_other_score, v.consensus_score);
    }
    REQUIRE(iso->consensus_score < min_other_score);
}

TEST_CASE("Hybrid 'robust' preset — thermal noise T~300K: pristine stays at 0 vacancies") {
    const int nc = 6;
    const double a = 2.87;
    Frame ref = make_bcc(nc, a);

    // σ ≈ 0.08 Å is a reasonable thermal RMS for Fe BCC at ~300 K.
    Frame dmg = ref;
    thermalize(dmg, 0.08, /*seed=*/1234);

    WignerSeitz ws;
    ws.build(ref, 1.6 * a);
    auto ws_res = ws.classify(dmg);
    // With moderate noise no atom should leave its cell; WS = 0 vacancies.
    REQUIRE(ws.vacancyCount() == 0);

    HybridParams p;
    p.applyPreset("robust");
    HybridVacancyDetector hvd(p);
    hvd.buildReference(ref, ws);
    hvd.detect(dmg, ws_res);
    REQUIRE(hvd.vacancyCount() == 0);
}
