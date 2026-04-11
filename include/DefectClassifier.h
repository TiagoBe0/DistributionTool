#pragma once
#include "AtomData.h"
#include "Statistics.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace DistTool {

/**
 * Reference descriptor set for a given material and temperature.
 *
 * Holds:
 *   - Mean DV  q̄(T)  from a thermalized perfect crystal.
 *   - Variance of distances (used in the Gaussian probability model, Eq. 6).
 *   - Optional reference DVs for known defect types (interstitial, vacancy-adj,
 *     type-a) used for secondary classification.
 */
struct ReferenceSet {
    std::vector<double> mean_dv;           // q̄(T)
    double              mean_dist = 0.0;   // <d>  of thermalized sample
    double              var_dist  = 0.0;   // <d²> — used in Eq. 6

    // Optional per-defect reference DVs (may be empty)
    std::vector<double> dv_interstitial;
    std::vector<double> dv_vacancy_adj;
    std::vector<double> dv_type_a;

    bool isSet() const { return !mean_dv.empty(); }
};

/**
 * DefectClassifier
 *
 * Pipeline:
 *   1. buildReference()  — compute q̄(T) from a pristine/thermalized frame.
 *   2. classify()        — label every atom in a damaged frame.
 *   3. findVacancies()   — locate vacant lattice sites via k-d brute-force.
 *
 * Classification rule (§2.3 of the paper):
 *   d^i = ||q̃^i − q̄(T)||
 *   If d^i < threshold  → Lattice
 *   Else compare d^i to reference DVs of known defects and assign nearest type.
 */
class DefectClassifier {
public:
    explicit DefectClassifier(double threshold = 0.15);

    // Build the reference set from a set of descriptor vectors
    // (all atoms in the pristine/thermalized frame).
    void buildReference(const std::vector<std::vector<double>>& ref_dvs);

    // Optionally add reference DVs for known defect types
    // to improve secondary classification accuracy.
    void setDefectReferences(
        const std::vector<double>& dv_interstitial,
        const std::vector<double>& dv_vacancy_adj,
        const std::vector<double>& dv_type_a = {});

    // Classify every atom in `frame`.  Populates dist_to_ref, defect_prob,
    // and defect_type on each atom.
    void classify(Frame& frame) const;

    // Sampling-grid vacancy detection (method of FaVaD, §2.3.2).
    // Places a uniform grid of points inside the simulation box.
    // Grid points whose nearest atom in `damaged` is farther than
    // dist_threshold [Å] are returned as vacancy/void positions.
    // Does NOT require a pristine reference frame.
    std::vector<std::array<double,3>> findVacanciesGrid(
        const Frame& damaged,
        double grid_spacing,       // grid point separation [Å]
        double dist_threshold) const;

    // Legacy: identify vacant sites from pristine lattice positions.
    // Kept for comparison; findVacanciesGrid is preferred.
    std::vector<std::array<double,3>> findVacancies(
        const Frame& pristine,
        const Frame& damaged,
        double dist_threshold) const;

    double           threshold()  const { return threshold_; }
    const ReferenceSet& reference() const { return ref_; }

private:
    double       threshold_;
    ReferenceSet ref_;

    DefectType secondaryClassify(const Atom& atom) const;
};

} // namespace DistTool
