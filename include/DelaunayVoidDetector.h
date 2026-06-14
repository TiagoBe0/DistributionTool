#pragma once
#include "AtomData.h"
#include <array>
#include <vector>

namespace DistTool {

/**
 * Reference-free vacancy detector based on periodic 3D Delaunay triangulation.
 *
 * A vacancy is a "hole" in the atomic point cloud — a region that should contain
 * an atom but doesn't.  Finding such holes requires no reference frame: a Delaunay
 * tetrahedron whose circumsphere is larger than the typical nearest-neighbour
 * distance sits in a void.  Adjacent void tetrahedra are clustered; each connected,
 * roughly-spherical cluster is counted as one point vacancy.
 *
 * Advantages over Wigner-Seitz:
 *   - No reference frame needed → works for any HEA composition or ordering.
 *   - Grain boundaries produce elongated clusters (high aspect ratio) and are
 *     filtered out automatically.
 *   - Uniaxial compression/expansion: the threshold scales with the local
 *     nearest-neighbour distance, not with a fixed reference lattice.
 *
 * Requires CGAL (libcgal-dev).  Built only when USE_CGAL is defined.
 */

struct DelaunayVoid {
    std::array<double,3> pos;   // cluster centroid (wrapped into simulation box)
    double circumradius;         // max circumsphere radius in cluster [Å]
    int    n_cells;              // Delaunay cells in this void cluster
    double extent;               // max pairwise distance between void-cell centres [Å]
    double aspect_ratio;         // sqrt(λ_max/λ_min): 1.0 = spherical, high = planar/linear
    double volume;               // summed volume of the void cells [Å³]
    int    n_vacancies;          // estimated vacancies in this cluster (≥1 if not extended)
    bool   is_extended;          // GB / dislocation / surface (extent ≫ d_nn) → not a vacancy
};

struct DelaunayVoidParams {
    // R_circumsphere / d_nn > threshold_ratio → cell flagged as void.
    // In a perfect BCC/FCC crystal the max circumratio is ~0.87; a monovacancy
    // creates cells with ratio ≥ 1.0.  Default 0.90 catches all vacancies with
    // very few false positives in the bulk.
    double threshold_ratio   = 0.90;

    // A void cluster whose linear extent exceeds extended_factor · d_nn is an
    // EXTENDED defect (grain-boundary sheet, dislocation-core tube, surface) and
    // is reported separately, not counted as a vacancy.  A mono/di/tri-vacancy
    // spans only a few d_nn, so this — NOT the aspect ratio — is the right
    // discriminator: a di-vacancy is "elongated" (high aspect ratio) yet tiny.
    double extended_factor   = 3.0;   // extent threshold = extended_factor · d_nn

    // Volume occupied by a single vacancy, used to split a compact void cluster
    // into a vacancy count (n_vac ≈ round(cluster_volume / vac_volume)).
    // -1 = auto: calibrated to the per-vacancy void volume of the lattice
    // (≈ vac_volume_factor · vol_per_atom; see .cpp).
    //
    // The factor 3.2 is the summed volume of the empty Delaunay cells around one
    // vacancy, in units of vol_per_atom — calibrated on a relaxed FeCrNi (BCC)
    // frame where Delaunay's per-cluster volumes were {mono ≈ 35 Å³, di ≈ 65 Å³}
    // against vol_per_atom ≈ 11 Å³, reproducing the WS count (6 mono + 2 di = 10)
    // exactly. May need re-calibration for FCC; override with --dv-vac-volume.
    double vac_volume        = -1.0;
    double vac_volume_factor = 3.2;   // void-cell volume per vacancy ≈ 3.2 · vol_per_atom

    // If > 0, use this value as d_nn instead of auto-detecting from frame density.
    double nn_override       = -1.0;
};

class DelaunayVoidDetector {
public:
    explicit DelaunayVoidDetector(const DelaunayVoidParams& p = {});

    // Calibrate d_nn from a pristine reference frame (optional).
    // When omitted, d_nn is estimated from the atom density of the damage frame.
    void calibrate(const Frame& ref_frame);

    // Detect void clusters in `frame` using periodic Delaunay triangulation.
    // Does NOT require a reference frame.
    std::vector<DelaunayVoid> detect(const Frame& frame) const;

    // Total estimated vacancies = Σ n_vacancies over non-extended clusters.
    int vacancyCount(const std::vector<DelaunayVoid>& voids) const;

    // Number of distinct vacancy clusters (non-extended voids).
    int clusterCount(const std::vector<DelaunayVoid>& voids) const;

    // Number of extended defects (GB / dislocation / surface voids).
    int extendedCount(const std::vector<DelaunayVoid>& voids) const;

    double nnRef() const { return nn_ref_; }

    // Estimate d_nn from atom number density (structure-independent, ~5% error).
    // Public so callers can preview the auto-detected value before running detect().
    static double estimateNN(const Frame& frame);

private:
    DelaunayVoidParams p_;
    double nn_ref_ = 0.0;  // 0 = auto-detect from frame density
};

} // namespace DistTool
