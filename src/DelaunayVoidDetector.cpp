#include "DelaunayVoidDetector.h"

// ── CGAL periodic Delaunay triangulation ─────────────────────────────────────
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

// ── Eigen (covariance / aspect ratio) ────────────────────────────────────────
#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace DistTool {
namespace {

// CGAL type aliases — kept in .cpp so CGAL headers never leak into the rest.
using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using GT  = CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>;
using PDT = CGAL::Periodic_3_Delaunay_triangulation_3<GT>;
using CGALPoint  = GT::Point_3;
using CellHandle = PDT::Cell_handle;

// Wrap coordinate into [lo, hi).
inline double wrapCoord(double v, double lo, double hi) {
    const double L = hi - lo;
    v = std::fmod(v - lo, L);
    if (v < 0.0) v += L;
    return v + lo;
}

// Wrap `p` relative to `ref`, choosing the periodic image closest to `ref`.
inline void wrapRelativeTo(double* p, const double* ref, const double* L) {
    for (int d = 0; d < 3; d++) {
        double diff = p[d] - ref[d];
        if      (diff >  0.5 * L[d]) p[d] -= L[d];
        else if (diff < -0.5 * L[d]) p[d] += L[d];
    }
}

// Integer key from a Cell_handle (raw pointer address).
inline std::size_t cellKey(CellHandle h) {
    return reinterpret_cast<std::size_t>(&*h);
}

// Analytic circumcenter of a tetrahedron (p0..p3).  Computed directly from the
// four vertices — which CGAL hands us in a single consistent spatial frame via
// periodic_point(c,i) — so the result never suffers the offset-wrapping bug that
// would arise from differencing canonical periodic points across a box boundary.
// Returns false for (near-)degenerate cells (slivers): their analytic
// circumcenter is numerically unstable, and a flat sliver is not a real void.
inline bool tetraCircum(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                        const Eigen::Vector3d& p2, const Eigen::Vector3d& p3,
                        Eigen::Vector3d& cc, double& r2) {
    const Eigen::Vector3d a = p1 - p0, b = p2 - p0, c = p3 - p0;
    const double denom = 2.0 * a.dot(b.cross(c));   // 2·(signed volume)·6
    if (std::abs(denom) < 1e-9) return false;
    const Eigen::Vector3d rel =
        ( a.squaredNorm() * b.cross(c)
        + b.squaredNorm() * c.cross(a)
        + c.squaredNorm() * a.cross(b) ) / denom;
    cc = p0 + rel;
    r2 = rel.squaredNorm();
    return true;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────

DelaunayVoidDetector::DelaunayVoidDetector(const DelaunayVoidParams& p)
    : p_(p), nn_ref_(p.nn_override > 0.0 ? p.nn_override : 0.0) {}

double DelaunayVoidDetector::estimateNN(const Frame& frame) {
    // Structure-independent estimate: d_nn ≈ 1.1 × (V/N)^(1/3).
    // Derivation: for both BCC and FCC the ratio d_nn / a^(1/3) ≈ 1.09-1.12
    // (BCC: d_nn = √3/2 · a, V/N = a³/2; FCC: d_nn = a/√2, V/N = a³/4).
    const SimBox& b = frame.box;
    const double vol = b.lx() * b.ly() * b.lz();
    const double vol_per_atom = vol / static_cast<double>(frame.atoms.size());
    return 1.1 * std::cbrt(vol_per_atom);
}

void DelaunayVoidDetector::calibrate(const Frame& ref_frame) {
    nn_ref_ = estimateNN(ref_frame);
}

std::vector<DelaunayVoid> DelaunayVoidDetector::detect(const Frame& frame) const {
    if (frame.atoms.empty())
        throw std::runtime_error("DelaunayVoidDetector::detect: empty frame");

    const SimBox& box = frame.box;
    const double xlo = box.xb[0], ylo = box.yb[0], zlo = box.zb[0];
    const double xhi = box.xb[1], yhi = box.yb[1], zhi = box.zb[1];
    const double L[3] = {box.lx(), box.ly(), box.lz()};
    const double lo[3] = {xlo, ylo, zlo};

    // Effective d_nn: from calibration, nn_override, or auto-estimate.
    const double d_nn = (nn_ref_ > 0.0)     ? nn_ref_
                      : (p_.nn_override > 0.0) ? p_.nn_override
                      : estimateNN(frame);

    const double r_thresh2 = std::pow(p_.threshold_ratio * d_nn, 2.0);

    // Optional phase timing (set DV_TIMING=1) to profile the triangulation cost.
    const bool dv_timing = std::getenv("DV_TIMING") != nullptr;
    using clk = std::chrono::steady_clock;
    auto t_mark = clk::now();
    auto lap = [&](const char* label) {
        if (!dv_timing) return;
        auto now = clk::now();
        const double ms = std::chrono::duration<double, std::milli>(now - t_mark).count();
        std::cerr << "[Delaunay timing] " << label << ": " << ms << " ms\n";
        t_mark = now;
    };

    // ── 1. Build periodic Delaunay triangulation ─────────────────────────────
    GT::Iso_cuboid_3 domain(xlo, ylo, zlo, xhi, yhi, zhi);
    PDT dt(domain);

    // Range insertion (NOT point-by-point): CGAL spatially sorts the whole batch
    // and inserts along a locality-preserving curve.  This is the difference
    // between seconds and hours for a 10⁶-atom frame — incremental one-at-a-time
    // insertion into a periodic triangulation is pathologically slow at this scale.
    std::vector<CGALPoint> pts;
    pts.reserve(frame.atoms.size());
    for (const auto& atom : frame.atoms) {
        // CGAL requires points strictly inside the domain.
        pts.emplace_back(wrapCoord(atom.x, xlo, xhi),
                         wrapCoord(atom.y, ylo, yhi),
                         wrapCoord(atom.z, zlo, zhi));
    }
    dt.insert(pts.begin(), pts.end(), /*is_large_point_set=*/true);
    lap("range insert");
    if (dv_timing)
        std::cerr << "[Delaunay timing]   1-cover=" << dt.is_1_cover()
                  << "  cells=" << dt.number_of_cells()
                  << "  vertices=" << dt.number_of_vertices() << "\n";

    // Periodic Delaunay starts in a 27-sheeted cover and collapses to a single
    // sheet only once that is provably a simplicial complex.  Iterating cells
    // while still in the 27-cover would count every tetrahedron 27×.  A full MD
    // frame (dense, box ≫ largest void) always reaches the 1-cover; bail loudly
    // otherwise rather than return a silently inflated count.
    if (!dt.is_1_cover())
        throw std::runtime_error(
            "DelaunayVoidDetector: periodic triangulation is in a 27-sheeted "
            "cover (simulation box too small or too few atoms). Reference-free "
            "detection requires a 1-sheeted cover.");

    // ── 2. Index all cells; flag those with circumsphere R > threshold ────────
    std::unordered_map<std::size_t, int> cell_idx;
    cell_idx.reserve(dt.number_of_cells());
    std::vector<CellHandle> all_cells;
    all_cells.reserve(dt.number_of_cells());

    for (auto cit = dt.cells_begin(); cit != dt.cells_end(); ++cit) {
        cell_idx[cellKey(cit)] = static_cast<int>(all_cells.size());
        all_cells.push_back(cit);
    }
    lap("index cells");

    const int n_cells = static_cast<int>(all_cells.size());
    std::vector<bool>   is_void(n_cells, false);
    std::vector<double> circ_r2(n_cells, 0.0);
    std::vector<double> cell_vol(n_cells, 0.0);
    std::vector<std::array<double,3>> circ_center(n_cells, {0,0,0});

    for (int i = 0; i < n_cells; i++) {
        CellHandle cit = all_cells[i];
        // Four vertices in this cell's consistent spatial frame (offsets applied).
        Eigen::Vector3d pv[4];
        for (int k = 0; k < 4; k++) {
            CGALPoint p = dt.construct_point(dt.periodic_point(cit, k));
            pv[k] = Eigen::Vector3d(CGAL::to_double(p.x()),
                                    CGAL::to_double(p.y()),
                                    CGAL::to_double(p.z()));
        }
        Eigen::Vector3d cc;
        double r2;
        if (!tetraCircum(pv[0], pv[1], pv[2], pv[3], cc, r2)) {
            is_void[i] = false;          // degenerate sliver → not a void
            continue;
        }
        circ_r2[i]     = r2;
        // Tetrahedron volume |det(p1-p0, p2-p0, p3-p0)| / 6.
        cell_vol[i]    = std::abs((pv[1]-pv[0]).dot((pv[2]-pv[0]).cross(pv[3]-pv[0]))) / 6.0;
        circ_center[i] = {wrapCoord(cc.x(), xlo, xhi),
                          wrapCoord(cc.y(), ylo, yhi),
                          wrapCoord(cc.z(), zlo, zhi)};
        is_void[i]     = (r2 > r_thresh2);
    }
    lap("circumcenters");

    // d_nn ⇒ per-vacancy void volume for splitting compact clusters into a count.
    const double vol_per_atom = std::pow(d_nn / 1.1, 3.0);   // inverse of estimateNN
    const double vac_vol = (p_.vac_volume > 0.0)
                         ? p_.vac_volume
                         : p_.vac_volume_factor * vol_per_atom;
    const double extent_thresh = p_.extended_factor * d_nn;

    // ── 3. BFS: connected components of void cells ───────────────────────────
    std::vector<bool> visited(n_cells, false);
    std::vector<DelaunayVoid> results;

    for (int seed = 0; seed < n_cells; seed++) {
        if (!is_void[seed] || visited[seed]) continue;

        std::queue<int> q;
        q.push(seed);
        visited[seed] = true;
        std::vector<int> cluster;

        while (!q.empty()) {
            int ci = q.front(); q.pop();
            cluster.push_back(ci);

            CellHandle cit = all_cells[ci];
            for (int f = 0; f < 4; f++) {
                CellHandle nbr = cit->neighbor(f);
                auto it = cell_idx.find(cellKey(nbr));
                if (it == cell_idx.end()) continue;
                int ni = it->second;
                if (!visited[ni] && is_void[ni]) {
                    visited[ni] = true;
                    q.push(ni);
                }
            }
        }

        // ── 4. Cluster properties ────────────────────────────────────────────
        const int nc = static_cast<int>(cluster.size());

        // Reuse the circumcenters/radii/volumes computed once in Step 2.
        std::vector<std::array<double,3>> centers;
        centers.reserve(nc);
        double max_r2 = 0.0, vol_sum = 0.0;

        for (int ci : cluster) {
            centers.push_back(circ_center[ci]);
            vol_sum += cell_vol[ci];
            if (circ_r2[ci] > max_r2) max_r2 = circ_r2[ci];
        }

        // PBC-aware unwrap: bring every centre into the same image as centre 0,
        // then take the centroid.  Unwrapped coords are reused for extent/shape.
        std::array<double,3> centroid = centers[0];
        for (int i = 1; i < nc; i++) {
            wrapRelativeTo(centers[i].data(), centers[0].data(), L);
            for (int d = 0; d < 3; d++)
                centroid[d] += centers[i][d];
        }
        for (int d = 0; d < 3; d++) centroid[d] /= nc;

        // Linear extent = max pairwise distance between void-cell centres [Å].
        // This is the physical size of the cluster and the discriminator between
        // a compact vacancy cluster and an extended defect — independent of the
        // aspect ratio, which is high even for a tiny di-vacancy.
        double extent = 0.0;
        for (int i = 0; i < nc; i++)
            for (int j = i + 1; j < nc; j++) {
                double d2 = 0.0;
                for (int d = 0; d < 3; d++) {
                    const double dd = centers[i][d] - centers[j][d];
                    d2 += dd * dd;
                }
                if (d2 > extent * extent) extent = std::sqrt(d2);
            }

        // Aspect ratio via covariance-matrix eigendecomposition (needs ≥ 4 points).
        double aspect_ratio = 1.0;
        if (nc >= 4) {
            Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
            for (auto& c : centers) {
                Eigen::Vector3d dp(c[0] - centroid[0],
                                   c[1] - centroid[1],
                                   c[2] - centroid[2]);
                cov += dp * dp.transpose();
            }
            cov /= nc;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov,
                                                                Eigen::EigenvaluesOnly);
            const auto& ev = eig.eigenvalues();  // ascending
            if (ev[0] > 1e-12)
                aspect_ratio = std::sqrt(ev[2] / ev[0]);
        }

        DelaunayVoid dv;
        for (int d = 0; d < 3; d++)
            dv.pos[d] = wrapCoord(centroid[d], lo[d], lo[d] + L[d]);
        dv.circumradius  = std::sqrt(max_r2);
        dv.n_cells       = nc;
        dv.extent        = extent;
        dv.aspect_ratio  = aspect_ratio;
        dv.volume        = vol_sum;
        dv.is_extended   = (extent > extent_thresh);
        // A compact cluster holds round(volume / per-vacancy void volume) vacancies,
        // floored at 1.  Extended defects are not assigned a vacancy count.
        dv.n_vacancies   = dv.is_extended
                         ? 0
                         : std::max(1, (int)std::lround(vol_sum / vac_vol));
        results.push_back(dv);
    }

    // Sort: vacancy clusters first, then by vacancy count descending.
    std::sort(results.begin(), results.end(),
              [](const DelaunayVoid& a, const DelaunayVoid& b) {
                  if (a.is_extended != b.is_extended)
                      return a.is_extended < b.is_extended;   // vacancies first
                  return a.n_vacancies > b.n_vacancies;
              });
    lap("BFS + cluster props");

    return results;
}

int DelaunayVoidDetector::vacancyCount(const std::vector<DelaunayVoid>& voids) const {
    int n = 0;
    for (const auto& v : voids) n += v.n_vacancies;   // 0 for extended defects
    return n;
}

int DelaunayVoidDetector::clusterCount(const std::vector<DelaunayVoid>& voids) const {
    int n = 0;
    for (const auto& v : voids) if (!v.is_extended) ++n;
    return n;
}

int DelaunayVoidDetector::extendedCount(const std::vector<DelaunayVoid>& voids) const {
    int n = 0;
    for (const auto& v : voids) if (v.is_extended) ++n;
    return n;
}

} // namespace DistTool
