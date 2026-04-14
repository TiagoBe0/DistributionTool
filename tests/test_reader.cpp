#include "test_harness.h"
#include "LammpsDumpReader.h"
#include <algorithm>
#include <string>

using namespace DistTool;

// TEST_DATA_DIR is injected by CMake (-DTEST_DATA_DIR=...)
#ifndef TEST_DATA_DIR
#  define TEST_DATA_DIR "."
#endif

static const std::string NANO_DUMP =
    std::string(TEST_DATA_DIR) + "/dump.nanoparticula.0";

TEST_CASE("LammpsDumpReader — parses nanoparticle dump: atom count") {
    LammpsDumpReader r(NANO_DUMP);
    REQUIRE(r.hasNext());
    Frame f = r.readNext();
    REQUIRE(f.size() == 30921);
}

TEST_CASE("LammpsDumpReader — parses nanoparticle dump: box bounds") {
    LammpsDumpReader r(NANO_DUMP);
    Frame f = r.readNext();
    // Box is [-75, 75]^3
    APPROX_EQ(f.box.xb[0], -75.0, 1e-6);
    APPROX_EQ(f.box.xb[1],  75.0, 1e-6);
    APPROX_EQ(f.box.yb[0], -75.0, 1e-6);
    APPROX_EQ(f.box.zb[1],  75.0, 1e-6);
}

TEST_CASE("LammpsDumpReader — parses nanoparticle dump: atom IDs are unique") {
    LammpsDumpReader r(NANO_DUMP);
    Frame f = r.readNext();
    // After sort-by-id, consecutive IDs must differ
    for (int i = 1; i < f.size(); ++i)
        REQUIRE(f.atoms[i].id != f.atoms[i-1].id);
}

TEST_CASE("LammpsDumpReader — parses nanoparticle dump: all positions inside box") {
    LammpsDumpReader r(NANO_DUMP);
    Frame f = r.readNext();
    const double xlo = f.box.xb[0], xhi = f.box.xb[1];
    const double ylo = f.box.yb[0], yhi = f.box.yb[1];
    const double zlo = f.box.zb[0], zhi = f.box.zb[1];
    for (const auto& a : f.atoms) {
        REQUIRE(a.x >= xlo && a.x <= xhi);
        REQUIRE(a.y >= ylo && a.y <= yhi);
        REQUIRE(a.z >= zlo && a.z <= zhi);
    }
}

TEST_CASE("LammpsDumpReader — hasNext is false after reading single-frame file") {
    LammpsDumpReader r(NANO_DUMP);
    r.readNext();
    REQUIRE(!r.hasNext());
}
