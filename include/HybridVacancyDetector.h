#pragma once
#include "AtomData.h"
#include "DefectClassifier.h"
#include "WignerSeitz.h"
#include <array>
#include <string>
#include <vector>

namespace DistTool {

/**
 * HybridVacancyDetector
 *
 * Super-method that combines several independent signals into a single
 * consensus score per candidate vacancy.  Designed to:
 *   - Reproduce Wigner-Seitz exactly when only the WS weight is enabled
 *     (default preset = "ws", reproduces WignerSeitz::vacancyCount()).
 *   - Reduce WS over-estimation in high-deformation regimes (PKA cascades)
 *     by penalising candidates where the "vacant" site has an atom in
 *     ballistic transit nearby.
 *
 * Signals (each squashed to [0,1] via sigmoid):
 *   1. WS occupancy (binary)                 — from WignerSeitz
 *   2. Soft-WS gaussian occupancy            — thermal-noise robust
 *   3. NN-excess (Voronoi proxy)             — neighbour-distance anomaly
 *   4. Density deficit (KDE)                 — reference-free void detector
 *   5. SOAP neighbor anomaly                 — uses dist_to_ref already populated
 *   6. Topology / coordination anomaly       — # neighbours vs ideal lattice
 *
 * Penalties (subtracted from the consensus):
 *   P1. Transit filter (ballistic atom near a "vacant" site is no defect)
 *   P2. Frenkel pair coherence (vacancy + nearby interstitial = single displaced atom)
 *
 * Consensus formula:
 *   score = (Σ_k w_k · s_k) / (Σ_k w_k)  −  w_transit · p_transit  −  w_frenkel · p_frenkel
 *
 * A candidate is reported as a vacancy iff `score >= accept_threshold` and
 * `p_transit < 0.5`.
 */

struct VacancySignal {
    double ws              = 0.0;
    double soft_ws         = 0.0;
    double voronoi_anomaly = 0.0;
    double density_deficit = 0.0;
    double soap_neighbor   = 0.0;
    double topology_anomaly= 0.0;
    double transit_penalty = 0.0;
    double frenkel_penalty = 0.0;
    double consensus_score = 0.0;
};

struct HybridVacancy {
    std::array<double,3> pos    {0,0,0};
    int    ref_site_idx         = -1;     // index into WignerSeitz::sites() (−1 if grid-only)
    double consensus_score      = 0.0;
    VacancySignal breakdown;
    bool   transit_filtered     = false;
    bool   in_recombination     = false;
};

struct HybridParams {
    // Per-signal weights.  Defaults reproduce WS exactly.
    double w_ws        = 1.0;
    double w_soft_ws   = 0.0;
    double w_voronoi   = 0.0;
    double w_density   = 0.0;
    double w_soap      = 0.0;
    double w_topology  = 0.0;
    double w_transit   = 0.0;
    double w_frenkel   = 0.0;

    double accept_threshold = 0.5;        // [0,1]

    // Physical parameters
    double thermal_sigma   = -1.0;        // -1 = auto (0.05 · nn_ref)
    double recomb_radius   = 3.3;         // Å, Fe-BCC default
    double transit_factor  = 1.5;         // < transit_factor · nn_ref ⇒ ballistic
    double grid_spacing    = 0.5;         // Å, for grid-vacancy seed candidates

    // Apply a named preset (overwrites all weights).
    //   "ws"        — w_ws=1, rest=0  (reproduces WignerSeitz exactly)
    //   "robust"    — multi-signal consensus + transit/Frenkel penalties
    //                 (for ballistic / peak-damage snapshots)
    //   "relaxed"   — multi-signal consensus, penalties OFF
    //                 (for cooled / relaxed final structures; reconciles with WS)
    //   "sensitive" — wider net, lower accept_threshold
    void applyPreset(const std::string& name);
};

class HybridVacancyDetector {
public:
    explicit HybridVacancyDetector(const HybridParams& p);

    /// Auto-fit nn_ref, vol_per_atom, coord_ref, σ_thermal from the pristine
    /// reference frame.  Stores const pointers — caller must keep them alive.
    void buildReference(const Frame& ref_frame,
                        const WignerSeitz& ws_built,
                        const ReferenceSet* soap_ref = nullptr);

    /// Run all signals over the damaged frame and return the consensus
    /// vacancy list (sorted by consensus_score, descending).
    std::vector<HybridVacancy> detect(
        const Frame& dmg_frame,
        const std::vector<WSAtomResult>& ws_results);

    /// Diagnostic accessors (after detect())
    const std::vector<HybridVacancy>& candidates() const { return cand_; }
    int vacancyCount()    const;   // accepted candidates
    int wsAgreeCount()    const;   // WS says vac AND hybrid accepts
    int wsOnlyCount()     const;   // WS says vac, hybrid filters out
    int hybridOnlyCount() const;   // hybrid accepts, WS does NOT say vac

    double nnRef()        const { return nn_ref_; }
    double sigmaThermal() const { return sigma_thermal_; }
    double coordRef()     const { return coord_ref_; }

private:
    HybridParams                p_;
    const Frame*                ref_           = nullptr;
    const WignerSeitz*          ws_            = nullptr;
    const ReferenceSet*         soap_ref_      = nullptr;
    double                      nn_ref_        = 0.0;
    double                      vol_per_atom_  = 0.0;
    double                      sigma_thermal_ = 0.0;
    double                      coord_ref_     = 8.0;   // BCC default
    std::vector<HybridVacancy>  cand_;
};

} // namespace DistTool
