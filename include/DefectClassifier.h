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
    double              var_dist  = 0.0;   // variance of distances

    // Chi-distribution parameters fitted from the reference distances (FaVaD Eq. 6).
    // k ≈ number of active DV components; σ is a scale factor.
    double              k_chi     = 2.0;
    double              sigma_chi = 1.0;

    // Optional per-defect reference DVs (may be empty)
    std::vector<double> dv_interstitial;
    std::vector<double> dv_vacancy_adj;
    std::vector<double> dv_type_a;

    bool isSet() const { return !mean_dv.empty(); }
};

/// One empty-space grid point found during vacancy detection.
struct VacancyPoint {
    std::array<double,3> pos;
    double d_near;   ///< distance to nearest atom [Å]
};

/// One physical vacancy (or void cluster) produced by merging nearby grid points.
struct VacancyCluster {
    std::array<double,3> center;     ///< position of the grid point with max d_near
    double               d_near_max; ///< peak void depth in this cluster [Å]
    int                  n_pts;      ///< number of grid points merged into this cluster
};

/**
 * DefectClassifier
 *
 * Pipeline:
 *   1. buildReference()        — compute q̄(T) from a pristine/thermalized frame.
 *   2. classify()              — label every atom in a damaged frame.
 *   3. findVacanciesGrid()     — locate vacant grid points (FaVaD §2.3.2).
 *   4. clusterVacancyPoints()  — merge grid points into individual vacancy events.
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
    // dist_threshold [Å] are returned as VacancyPoint objects that include
    // the actual d_near value.  Does NOT require a pristine reference frame.
    std::vector<VacancyPoint> findVacanciesGrid(
        const Frame& damaged,
        double grid_spacing,       // grid point separation [Å]
        double dist_threshold) const;

    // Greedy clustering of vacancy grid points into individual vacancy events
    // (FaVaD §2.3.2 iterative algorithm).
    // Seeds clusters from the grid point with the largest d_near, absorbs all
    // unassigned points within r_cluster [Å], then repeats until done.
    // box is used to apply minimum-image PBC when comparing grid point distances;
    // pass box.periodic = {false,false,false} to disable PBC (e.g. nanoparticles).
    // Returns one VacancyCluster per physical vacancy/void detected.
    std::vector<VacancyCluster> clusterVacancyPoints(
        const std::vector<VacancyPoint>& pts,
        double r_cluster,
        const SimBox& box) const;

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
