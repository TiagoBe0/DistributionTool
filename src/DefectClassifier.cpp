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

        // Defect probability  1 − P(d^i | T)
        double p_lattice = Statistics::latticeProbability(
            atom.dist_to_ref, ref_.mean_dist, ref_.var_dist);
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

// ── Vacancy identification (Appendix A of the paper) ─────────────────────────
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
