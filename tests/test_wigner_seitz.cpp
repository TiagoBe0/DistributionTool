#include "test_harness.h"
#include "AtomData.h"
#include "WignerSeitz.h"
#include <cmath>
#include <vector>

using namespace DistTool;

// ── Helpers ───────────────────────────────────────────────────────────────────

// Build a non-periodic box big enough to hold all atoms.
static SimBox make_box(double lx, double ly, double lz) {
    SimBox b;
    b.xb = {0, lx}; b.yb = {0, ly}; b.zb = {0, lz};
    b.periodic[0] = b.periodic[1] = b.periodic[2] = false;
    return b;
}

// Chain of N atoms along x at positions (i*a + a/2, 5, 5).
static Frame make_chain(int N, double a) {
    Frame f;
    f.box = make_box(N * a + 1.0, 10.0, 10.0);
    for (int i = 0; i < N; ++i) {
        Atom at;
        at.id = i + 1; at.type = 1;
        at.x = i * a + a * 0.5; at.y = 5.0; at.z = 5.0;
        f.atoms.push_back(at);
    }
    return f;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("WignerSeitz::build — r_cut <= 0 throws") {
    WignerSeitz ws;
    Frame f = make_chain(4, 3.0);
    CHECK_THROWS(ws.build(f, 0.0));
    CHECK_THROWS(ws.build(f, -1.0));
}

TEST_CASE("WignerSeitz::build — sites count matches reference frame") {
    WignerSeitz ws;
    Frame ref = make_chain(10, 3.0);
    ws.build(ref, 4.0);
    REQUIRE(static_cast<int>(ws.sites().size()) == 10);
}

TEST_CASE("WignerSeitz — perfect crystal: all occupancies 1, all Lattice") {
    const int N = 6;
    const double a = 3.0;
    WignerSeitz ws;
    Frame ref = make_chain(N, a);
    ws.build(ref, a * 1.5);

    Frame dmg = ref;  // identical positions
    auto results = ws.classify(dmg);

    REQUIRE(static_cast<int>(results.size()) == N);
    for (int i = 0; i < N; ++i) {
        REQUIRE(results[i].ref_site_idx >= 0);
        REQUIRE(results[i].ws_occ == 1);
        REQUIRE(!results[i].is_interstitial);
        APPROX_EQ(results[i].dist, 0.0, 1e-10);
    }
    REQUIRE(ws.vacancyCount() == 0);
    REQUIRE(ws.interstitialCount() == 0);
}

TEST_CASE("WignerSeitz — one vacancy: site occupancy 0, correct counts") {
    const int N = 5;
    const double a = 3.0;
    WignerSeitz ws;
    Frame ref = make_chain(N, a);
    ws.build(ref, a * 1.5);

    // Remove middle atom (index 2 = site at 2*a + a/2 = 7.5)
    Frame dmg = ref;
    dmg.atoms.erase(dmg.atoms.begin() + 2);
    REQUIRE(dmg.atoms.size() == static_cast<size_t>(N - 1));

    auto results = ws.classify(dmg);
    REQUIRE(static_cast<int>(results.size()) == N - 1);

    // All remaining atoms should map to their own sites (occ=1)
    for (int i = 0; i < N - 1; ++i) {
        REQUIRE(results[i].ws_occ == 1);
        REQUIRE(!results[i].is_interstitial);
    }

    REQUIRE(ws.vacancyCount() == 1);
    REQUIRE(ws.interstitialCount() == 0);

    // The vacant site should be site index 2 (x = 7.5)
    const auto& sites = ws.sites();
    APPROX_EQ(sites[2].x, 7.5, 1e-10);
    REQUIRE(sites[2].occupancy == 0);
}

TEST_CASE("WignerSeitz — one interstitial: crowded cell, correct counts") {
    const int N = 5;
    const double a = 3.0;
    WignerSeitz ws;
    Frame ref = make_chain(N, a);
    ws.build(ref, a * 1.5);

    // Add extra atom close to site 0 (x=1.5) → maps to same site
    Frame dmg = ref;
    Atom extra; extra.id = 99; extra.type = 1;
    extra.x = 1.5 + 0.2; extra.y = 5.0; extra.z = 5.0;
    dmg.atoms.push_back(extra);

    auto results = ws.classify(dmg);
    REQUIRE(static_cast<int>(results.size()) == N + 1);

    // Two atoms should map to site 0 → occ=2, both is_interstitial=true
    int crowded_count = 0;
    for (const auto& r : results)
        if (r.ref_site_idx == 0) ++crowded_count;
    REQUIRE(crowded_count == 2);

    for (const auto& r : results) {
        if (r.ref_site_idx == 0) {
            REQUIRE(r.ws_occ == 2);
            REQUIRE(r.is_interstitial);
        } else {
            REQUIRE(r.ws_occ == 1);
            REQUIRE(!r.is_interstitial);
        }
    }

    REQUIRE(ws.vacancyCount() == 0);
    REQUIRE(ws.interstitialCount() == 1);
}

TEST_CASE("WignerSeitz — Frenkel pair: 1 vacancy + 1 interstitial") {
    // Atom displaced from site 2 (x=7.5) to near site 0 (x=1.5) → site 2 vacant,
    // site 0 crowded (occ=2).
    const int N = 5;
    const double a = 3.0;
    WignerSeitz ws;
    Frame ref = make_chain(N, a);
    ws.build(ref, a * 1.5);

    Frame dmg = ref;
    // Move atom at index 2 (x=7.5) to (1.7, 5, 5) — closer to site 0 (1.5) than site 1 (4.5)
    dmg.atoms[2].x = 1.7;

    auto results = ws.classify(dmg);

    REQUIRE(ws.vacancyCount() == 1);
    REQUIRE(ws.interstitialCount() == 1);

    // Site 2 (x=7.5) is vacant
    REQUIRE(ws.sites()[2].occupancy == 0);
    // Site 0 (x=1.5) is doubly occupied
    REQUIRE(ws.sites()[0].occupancy == 2);
}

TEST_CASE("WignerSeitz — classify is idempotent (call twice gives same result)") {
    const int N = 6;
    WignerSeitz ws;
    Frame ref = make_chain(N, 3.0);
    ws.build(ref, 5.0);

    Frame dmg = ref;
    dmg.atoms[0].x += 0.1;  // small displacement

    auto r1 = ws.classify(dmg);
    auto r2 = ws.classify(dmg);

    REQUIRE(r1.size() == r2.size());
    for (int i = 0; i < N; ++i) {
        REQUIRE(r1[i].ref_site_idx == r2[i].ref_site_idx);
        REQUIRE(r1[i].ws_occ       == r2[i].ws_occ);
        APPROX_EQ(r1[i].dist,         r2[i].dist, 1e-12);
    }
}

TEST_CASE("WignerSeitz — sites vacancyCount=0 interstitialCount=0 before classify") {
    // After build() but before classify(), occupancy field is 0 for all sites.
    // vacancyCount() counts those as vacancies, which is N — this is a known
    // "uncalled classify" state.  Test that classify() resets it properly.
    const int N = 4;
    WignerSeitz ws;
    Frame ref = make_chain(N, 3.0);
    ws.build(ref, 5.0);

    // After classify with identical frame: all occ=1
    auto results = ws.classify(ref);
    REQUIRE(ws.vacancyCount()      == 0);
    REQUIRE(ws.interstitialCount() == 0);
    for (const auto& r : results)
        REQUIRE(r.ws_occ == 1);
}

TEST_CASE("WignerSeitz — nearest site with atom displaced beyond 1NN triggers fallback") {
    // 3 atoms in a line at x=1.5, 4.5, 7.5 (a=3.0). Atom at x=4.5 is moved to
    // x=9.5 (> r_cut=3.5). The fallback must find site index 2 (x=7.5) as nearest.
    const int N = 3;
    const double a = 3.0;
    WignerSeitz ws;
    Frame ref = make_chain(N, a);
    ws.build(ref, a * 1.1);  // r_cut=3.3: tighter than 9.5-7.5=2.0 but triggers fallback

    Frame dmg = ref;
    dmg.atoms[1].x = 9.5;  // moved well outside the 27-cell shell from its cell

    auto results = ws.classify(dmg);

    // Atom at 9.5 should map to site 2 (x=7.5, distance=2.0) not site 1 (x=4.5, dist=5.0)
    REQUIRE(results[1].ref_site_idx == 2);
    APPROX_EQ(results[1].dist, 2.0, 1e-10);

    // Site 1 (x=4.5) has no occupant → vacancy
    REQUIRE(ws.sites()[1].occupancy == 0);
    // Site 2 (x=7.5) has 2 occupants
    REQUIRE(ws.sites()[2].occupancy == 2);

    REQUIRE(ws.vacancyCount()      == 1);
    REQUIRE(ws.interstitialCount() == 1);
}
