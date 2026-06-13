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
    std::array<double,3> pos;   // vacancy centroid (wrapped into simulation box)
    double circumradius;         // max circumsphere radius in cluster [Å]
    int    n_cells;              // Delaunay cells in this void cluster
    double aspect_ratio;         // sqrt(λ_max/λ_min): 1.0 = spherical, high = planar
    bool   is_point_defect;     // passes size + shape filters → likely a point vacancy
};

struct DelaunayVoidParams {
    // R_circumsphere / d_nn > threshold_ratio → cell flagged as void.
    // In a perfect BCC/FCC crystal the max circumratio is ~0.87; a monovacancy
    // creates cells with ratio ≥ 1.0.  Default 0.90 catches all vacancies with
    // very few false positives in the bulk.
    double threshold_ratio   = 0.90;

    // Point vacancy clusters span at most this many Delaunay cells.
    // Extended defects (GB voids, dislocation loops) span far more.
    int    max_cluster_cells = 50;

    // Point vacancies have aspect ratio below this.
    // Planar GB voids have aspect_ratio >> 5.
    double max_aspect_ratio  = 5.0;

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

    // Count clusters that pass the point-defect filters.
    int vacancyCount(const std::vector<DelaunayVoid>& voids) const;

    double nnRef() const { return nn_ref_; }

    // Estimate d_nn from atom number density (structure-independent, ~5% error).
    // Public so callers can preview the auto-detected value before running detect().
    static double estimateNN(const Frame& frame);

private:
    DelaunayVoidParams p_;
    double nn_ref_ = 0.0;  // 0 = auto-detect from frame density
};

} // namespace DistTool
