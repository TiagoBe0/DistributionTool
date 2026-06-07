#pragma once
#include "AtomData.h"
#include "Statistics.h"
#include <map>
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

    // Reference distances sorted ascending — basis for the non-parametric
    // empirical defect probability and for percentile-derived thresholds.
    std::vector<double> sorted_dists;

    // Primary classification threshold for this set (-1 = use the global
    // classifier threshold). Set when a percentile threshold is requested.
    double              threshold = -1.0;

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

    // Enable per-species references: a separate q̄(T), chi-fit and distance
    // baseline is built and used for each atom `type`. Essential for chemically
    // disordered systems (e.g. high-entropy alloys), where a single global
    // reference mixes inequivalent local environments and collapses the chi-fit.
    // Must be called BEFORE buildReference(). Default: off (single global ref).
    void setPerSpecies(bool enabled) { per_species_ = enabled; }
    bool perSpecies() const { return per_species_; }

    // Use the non-parametric empirical CDF (percentile rank) for defect_prob
    // instead of the parametric chi model. Robust to non-chi distance
    // distributions (HEAs / disordered alloys). Default: off (chi model).
    void setEmpiricalProb(bool enabled) { empirical_prob_ = enabled; }
    bool empiricalProb() const { return empirical_prob_; }

    // Derive the primary classification threshold from the p-th percentile
    // (p ∈ [0,1]) of the reference distance distribution — per species when
    // per-species mode is on. Negative value (default) keeps the fixed
    // user-supplied threshold. Must be set BEFORE buildReference().
    void setThresholdPercentile(double p) { threshold_pct_ = p; }
    double thresholdPercentile() const { return threshold_pct_; }

    // Build the reference set from atoms in the pristine/thermalized frame.
    // DVs are read directly from atom.dv — no intermediate copy is made.
    // When per-species is enabled, one ReferenceSet is built per atom type
    // (in addition to the global reference, kept as a fallback/summary).
    void buildReference(const std::vector<Atom>& atoms);

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

    double           threshold()  const { return threshold_; }
    const ReferenceSet& reference() const { return ref_; }

    // Per-species references (empty unless per-species mode was enabled before
    // buildReference()). Keyed by atom type.
    const std::map<int, ReferenceSet>& referencesByType() const { return ref_by_type_; }

private:
    double       threshold_;
    bool         per_species_   = false;
    bool         empirical_prob_ = false;
    double       threshold_pct_  = -1.0;
    ReferenceSet ref_;                        // global reference / fallback
    std::map<int, ReferenceSet> ref_by_type_; // per-type references

    DefectType secondaryClassify(const Atom& atom) const;
};

} // namespace DistTool
