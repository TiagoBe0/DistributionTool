#pragma once
#include <array>
#include <string>
#include <vector>

namespace DistTool {

// ── Defect classification labels ─────────────────────────────────────────────
enum class DefectType : int {
    Lattice        = 0,  // Regular lattice site (undistorted)
    Interstitial   = 1,  // Self-interstitial atom
    VacancyAdjacent = 2, // Atom next to a single vacancy
    TypeA          = 3,  // Novel/unexpected defect (found via PCA)
    Unknown        = 4   // Distorted but unclassified
};

inline const char* defectName(DefectType t) {
    switch (t) {
        case DefectType::Lattice:          return "Lattice";
        case DefectType::Interstitial:     return "Interstitial";
        case DefectType::VacancyAdjacent:  return "VacancyAdj";
        case DefectType::TypeA:            return "TypeA";
        default:                           return "Unknown";
    }
}

// ── Per-atom data ─────────────────────────────────────────────────────────────
struct Atom {
    int    id   = 0;
    int    type = 1;
    double x = 0, y = 0, z = 0;

    // Descriptor vector q̃^i  (SOAP power spectrum, normalised)
    std::vector<double> dv;

    // Distance to mean reference DV:  d^i = ||q̃^i - q̄(T)||
    double dist_to_ref = 0.0;

    // 1 − P(q̃^i | T)  — probability of being in a distorted environment
    double defect_prob = 0.0;

    DefectType defect_type = DefectType::Lattice;
};

// ── Simulation box ────────────────────────────────────────────────────────────
struct SimBox {
    std::array<double, 2> xb{0, 1};  // [xlo, xhi]
    std::array<double, 2> yb{0, 1};
    std::array<double, 2> zb{0, 1};
    bool periodic[3] = {true, true, true};

    double lx() const { return xb[1] - xb[0]; }
    double ly() const { return yb[1] - yb[0]; }
    double lz() const { return zb[1] - zb[0]; }
};

// ── One snapshot from an MD trajectory ───────────────────────────────────────
struct Frame {
    int              timestep = 0;
    SimBox           box;
    std::vector<Atom> atoms;

    int size() const { return static_cast<int>(atoms.size()); }
};

} // namespace DistTool
