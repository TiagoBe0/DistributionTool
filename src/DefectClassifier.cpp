#include "DefectClassifier.h"
#include "CellList.h"
#include <algorithm>
#include <cmath>
#include <limits>
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
    Statistics::distanceStats(ref_dvs, ref_.mean_dv,
                               ref_.mean_dist, ref_.var_dist);

    // Fit chi-distribution parameters from the reference distance distribution.
    std::vector<double> dists;
    dists.reserve(ref_dvs.size());
    for (const auto& dv : ref_dvs)
        dists.push_back(Statistics::euclidean(dv, ref_.mean_dv));
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
// Points farther than dist_threshold from every atom are collected as
// vacancy/void positions.
std::vector<std::array<double,3>> DefectClassifier::findVacanciesGrid(
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

    std::vector<std::array<double,3>> vacancies;

    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
    for (int iz = 0; iz < nz; ++iz) {
        const double x = box.xb[0] + (ix + 0.5) * dx;
        const double y = box.yb[0] + (iy + 0.5) * dy;
        const double z = box.zb[0] + (iz + 0.5) * dz;
        if (cl.nearestDist2FromPoint(x, y, z) > dt2)
            vacancies.push_back({x, y, z});
    }

    return vacancies;
}

// ── Legacy: vacancy identification from pristine lattice positions ────────────
// Builds a cell list from the damaged frame for O(N) nearest-site lookup.
// Returns pristine lattice positions whose nearest atom in `damaged` exceeds
// `dist_threshold`.
std::vector<std::array<double,3>> DefectClassifier::findVacancies(
    const Frame& pristine,
    const Frame& damaged,
    double dist_threshold) const
{
    // Build a cell list from the damaged frame — O(N) construction,
    // O(27·k) per query where k = average atoms per cell.
    CellList cl;
    cl.build(damaged, dist_threshold);
    const double dt2 = dist_threshold * dist_threshold;

    std::vector<std::array<double,3>> vacancies;
    vacancies.reserve(pristine.atoms.size() / 16); // rough pre-alloc

    for (const auto& ref_atom : pristine.atoms) {
        if (cl.nearestDist2FromPoint(ref_atom.x, ref_atom.y, ref_atom.z) > dt2)
            vacancies.push_back({ref_atom.x, ref_atom.y, ref_atom.z});
    }

    return vacancies;
}

} // namespace DistTool
