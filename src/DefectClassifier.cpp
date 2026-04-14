#include "DefectClassifier.h"
#include "CellList.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace DistTool {

DefectClassifier::DefectClassifier(double threshold)
    : threshold_(threshold)
{}

// ── Build reference from pristine DVs ────────────────────────────────────────
void DefectClassifier::buildReference(
    const std::vector<std::vector<double>>& ref_dvs)
{
    if (ref_dvs.empty())
        throw std::invalid_argument("buildReference: empty DV set");

    ref_.mean_dv = Statistics::mean(ref_dvs);

    // Compute distances once; reuse for mean/var and chi-fit.
    std::vector<double> dists;
    dists.reserve(ref_dvs.size());
    for (const auto& dv : ref_dvs)
        dists.push_back(Statistics::euclidean(dv, ref_.mean_dv));

    const int nd = static_cast<int>(dists.size());
    ref_.mean_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / nd;
    double var = 0.0;
    for (double d : dists) { const double delta = d - ref_.mean_dist; var += delta * delta; }
    ref_.var_dist = var / nd;

    Statistics::fitChiParams(dists, ref_.k_chi, ref_.sigma_chi);
}

void DefectClassifier::setDefectReferences(
    const std::vector<double>& dv_interstitial,
    const std::vector<double>& dv_vacancy_adj,
    const std::vector<double>& dv_type_a)
{
    ref_.dv_interstitial = dv_interstitial;
    ref_.dv_vacancy_adj  = dv_vacancy_adj;
    ref_.dv_type_a       = dv_type_a;
}

// ── Secondary classification by nearest reference DV ─────────────────────────
DefectType DefectClassifier::secondaryClassify(const Atom& atom) const {
    auto dist = [&](const std::vector<double>& ref) -> double {
        if (ref.empty()) return std::numeric_limits<double>::max();
        return Statistics::euclidean(atom.dv, ref);
    };

    double di = dist(ref_.dv_interstitial);
    double dv = dist(ref_.dv_vacancy_adj);
    double da = dist(ref_.dv_type_a);

    if (di == std::numeric_limits<double>::max() &&
        dv == std::numeric_limits<double>::max() &&
        da == std::numeric_limits<double>::max())
        return DefectType::Unknown;

    if (di <= dv && di <= da) return DefectType::Interstitial;
    if (dv <= da)              return DefectType::VacancyAdjacent;
    return DefectType::TypeA;
}

// ── Classify all atoms in a damaged frame ────────────────────────────────────
void DefectClassifier::classify(Frame& frame) const {
    if (!ref_.isSet())
        throw std::runtime_error("classify: reference not built yet");

    for (auto& atom : frame.atoms) {
        // Primary distance metric (Eq. 5)
        atom.dist_to_ref = Statistics::euclidean(atom.dv, ref_.mean_dv);

        // Defect probability  1 − P(d^i | k, σ)  using chi-distribution model
        double p_lattice = Statistics::chiProbability(
            atom.dist_to_ref, ref_.k_chi, ref_.sigma_chi);
        atom.defect_prob = 1.0 - p_lattice;

        // Primary decision
        if (atom.dist_to_ref < threshold_) {
            atom.defect_type = DefectType::Lattice;
        } else {
            // Secondary: nearest known defect reference
            atom.defect_type = secondaryClassify(atom);
        }
    }
}

// ── Sampling-grid vacancy detection (FaVaD §2.3.2) ───────────────────────────
// Creates a Nx × Ny × Nz uniform grid within the simulation box.
// Each grid point is queried against the damaged frame's cell list.
// Points farther than dist_threshold from every atom are returned as
// VacancyPoint objects that include the actual d_near value (not just a flag).
std::vector<VacancyPoint> DefectClassifier::findVacanciesGrid(
    const Frame& damaged,
    double grid_spacing,
    double dist_threshold) const
{
    if (grid_spacing   <= 0.0) throw std::invalid_argument("grid_spacing must be > 0");
    if (dist_threshold <= 0.0) throw std::invalid_argument("dist_threshold must be > 0");

    const SimBox& box = damaged.box;
    const double Lx = box.lx(), Ly = box.ly(), Lz = box.lz();

    const int nx = std::max(1, static_cast<int>(std::ceil(Lx / grid_spacing)));
    const int ny = std::max(1, static_cast<int>(std::ceil(Ly / grid_spacing)));
    const int nz = std::max(1, static_cast<int>(std::ceil(Lz / grid_spacing)));

    const double dx = Lx / nx;
    const double dy = Ly / ny;
    const double dz = Lz / nz;
    const double dt2 = dist_threshold * dist_threshold;

    CellList cl;
    cl.build(damaged, dist_threshold);

    std::vector<VacancyPoint> vacancies;

    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
    for (int iz = 0; iz < nz; ++iz) {
        const double x = box.xb[0] + (ix + 0.5) * dx;
        const double y = box.yb[0] + (iy + 0.5) * dy;
        const double z = box.zb[0] + (iz + 0.5) * dz;
        const double d2 = cl.nearestDist2FromPoint(x, y, z);
        if (d2 > dt2)
            vacancies.push_back({{x, y, z}, std::sqrt(d2)});
    }

    return vacancies;
}

// ── Greedy vacancy clustering (FaVaD §2.3.2 iterative algorithm) ─────────────
// Seeds each cluster from the unassigned grid point with the largest d_near
// (deepest void), absorbs all unassigned points within r_cluster, repeats.
// Applies minimum-image PBC using the supplied SimBox.
// Returns one VacancyCluster per physical vacancy/void event.
//
// Complexity: O(M log M) sort + O(M) spatial grid build + O(M · 27k) clustering,
// where k = average vacancy points per grid cell.  Was O(M²) before.
std::vector<VacancyCluster> DefectClassifier::clusterVacancyPoints(
    const std::vector<VacancyPoint>& pts,
    double r_cluster,
    const SimBox& box) const
{
    if (r_cluster <= 0.0) throw std::invalid_argument("r_cluster must be > 0");

    const int n = static_cast<int>(pts.size());
    if (n == 0) return {};

    const double Lx = box.lx(), Ly = box.ly(), Lz = box.lz();

    // ── Spatial grid over vacancy points (same idea as CellList) ─────────────
    // Each cell side = r_cluster so the 3×3×3 shell covers the full search radius.
    const int gnx = std::max(1, static_cast<int>(std::floor(Lx / r_cluster)));
    const int gny = std::max(1, static_cast<int>(std::floor(Ly / r_cluster)));
    const int gnz = std::max(1, static_cast<int>(std::floor(Lz / r_cluster)));
    const double gcx = Lx / gnx, gcy = Ly / gny, gcz = Lz / gnz;
    const double xlo = box.xb[0], ylo = box.yb[0], zlo = box.zb[0];

    auto cellOf = [&](const std::array<double,3>& pos, int& ix, int& iy, int& iz) {
        ix = std::min(gnx-1, std::max(0, static_cast<int>((pos[0] - xlo) / gcx)));
        iy = std::min(gny-1, std::max(0, static_cast<int>((pos[1] - ylo) / gcy)));
        iz = std::min(gnz-1, std::max(0, static_cast<int>((pos[2] - zlo) / gcz)));
    };

    std::vector<std::vector<int>> grid(gnx * gny * gnz);
    for (int i = 0; i < n; ++i) {
        int ix, iy, iz;
        cellOf(pts[i].pos, ix, iy, iz);
        grid[ix * gny * gnz + iy * gnz + iz].push_back(i);
    }

    // Minimum-image displacement along one axis.
    auto mi = [](double d, double L, bool periodic) -> double {
        if (!periodic || L <= 0.0) return d;
        return d - L * std::round(d / L);
    };

    // Sort indices by d_near descending so we always seed from the deepest void.
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b){
        return pts[a].d_near > pts[b].d_near;
    });

    const double rc2 = r_cluster * r_cluster;
    std::vector<bool> assigned(n, false);
    std::vector<VacancyCluster> clusters;

    for (int ii = 0; ii < n; ++ii) {
        const int i = order[ii];
        if (assigned[i]) continue;

        // Seed new cluster from the deepest remaining point.
        VacancyCluster cl;
        cl.center     = pts[i].pos;
        cl.d_near_max = pts[i].d_near;
        cl.n_pts      = 1;
        assigned[i]   = true;

        // Query the 3×3×3 cell neighbourhood of the seed.
        int six, siy, siz;
        cellOf(pts[i].pos, six, siy, siz);

        for (int dix = -1; dix <= 1; ++dix)
        for (int diy = -1; diy <= 1; ++diy)
        for (int diz = -1; diz <= 1; ++diz) {
            int jx = six + dix, jy = siy + diy, jz = siz + diz;
            // Periodic: wrap; non-periodic: skip out-of-range shells.
            if (box.periodic[0]) jx = ((jx % gnx) + gnx) % gnx;
            else if (jx < 0 || jx >= gnx) continue;
            if (box.periodic[1]) jy = ((jy % gny) + gny) % gny;
            else if (jy < 0 || jy >= gny) continue;
            if (box.periodic[2]) jz = ((jz % gnz) + gnz) % gnz;
            else if (jz < 0 || jz >= gnz) continue;

            for (int j : grid[jx * gny * gnz + jy * gnz + jz]) {
                if (assigned[j]) continue;
                double ddx = mi(pts[j].pos[0] - cl.center[0], Lx, box.periodic[0]);
                double ddy = mi(pts[j].pos[1] - cl.center[1], Ly, box.periodic[1]);
                double ddz = mi(pts[j].pos[2] - cl.center[2], Lz, box.periodic[2]);
                if (ddx*ddx + ddy*ddy + ddz*ddz <= rc2) {
                    assigned[j] = true;
                    ++cl.n_pts;
                }
            }
        }

        clusters.push_back(cl);
    }

    return clusters;
}

} // namespace DistTool
