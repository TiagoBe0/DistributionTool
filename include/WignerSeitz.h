#pragma once
#include "AtomData.h"
#include <vector>

namespace DistTool {

/// One site from the reference (pristine) frame.
struct WSSite {
    int    ref_atom_id;          ///< ID of the reference atom
    double x, y, z;              ///< Reference lattice position [Å]
    int    occupancy = 0;        ///< Damaged atoms assigned to this WS cell
    double min_dist  = 0.0;      ///< Distance from nearest occupant to this site [Å]
};

/// Per damaged-atom result from WS classification.
struct WSAtomResult {
    int    ref_site_idx  = -1;   ///< Index into WignerSeitz::sites() (-1 if unassigned)
    double dist          = 0.0;  ///< Distance to the assigned reference site [Å]
    int    ws_occ        = 0;    ///< Occupancy of the assigned WS cell (1=normal, ≥2=crowded)
    bool   is_interstitial = false; ///< True when ws_occ ≥ 2 (extra atom in cell)
};

/**
 * WignerSeitz
 *
 * Classifies atoms in a damaged frame by assigning each one to the nearest
 * site in a reference (pristine) frame — the standard Wigner-Seitz method.
 *
 * Usage:
 *   WignerSeitz ws;
 *   ws.build(ref_frame, r_cut);          // build spatial index from reference
 *   auto results = ws.classify(dmg_frame); // assign and count occupancies
 *
 * Defect counting:
 *   Vacancies      = reference sites with occupancy == 0
 *   Interstitials  = sum of (occupancy - 1) over all sites with occupancy >= 2
 *
 * Per-atom label:
 *   Lattice        = assigned WS cell has occupancy == 1
 *   Interstitial   = assigned WS cell has occupancy >= 2
 *
 * Complexity: O(N_ref) build, O(N_dmg * 27k) classify.
 */
class WignerSeitz {
public:
    /// Build the spatial index from a pristine reference frame.
    /// r_cut must be larger than the maximum first-neighbour distance so that
    /// the 3×3×3 cell shell always covers the nearest reference site.
    void build(const Frame& ref_frame, double r_cut);

    /// Assign every atom in dmg_frame to its nearest reference site.
    /// Computes WS occupancies and returns one WSAtomResult per atom
    /// in the same order as dmg_frame.atoms.
    std::vector<WSAtomResult> classify(const Frame& dmg_frame);

    /// Reference sites with occupancy filled in after classify().
    const std::vector<WSSite>& sites() const { return sites_; }

    int vacancyCount() const;       ///< Sites with occupancy == 0
    int interstitialCount() const;  ///< Sum of max(0, occ-1) over all sites

private:
    std::vector<WSSite> sites_;

    // Spatial grid for O(N) nearest-site lookup
    SimBox ref_box_;
    double r_cut_ = 0.0;
    double Lx_ = 0, Ly_ = 0, Lz_ = 0;
    double xlo_ = 0, ylo_ = 0, zlo_ = 0;
    int    nx_ = 1, ny_ = 1, nz_ = 1;
    double cx_ = 1, cy_ = 1, cz_ = 1;
    std::vector<std::vector<int>> cells_;

    int  nearestSite(double qx, double qy, double qz, double& dist_out) const;
    void posCell(double x, double y, double z, int& ix, int& iy, int& iz) const;
    int  cellIdx(int ix, int iy, int iz) const { return ix * ny_ * nz_ + iy * nz_ + iz; }
};

} // namespace DistTool
