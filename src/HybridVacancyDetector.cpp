#include "HybridVacancyDetector.h"
#include "CellList.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace DistTool {

// ─────────────────────────────────────────────────────────────────────────────
// Local helpers
// ─────────────────────────────────────────────────────────────────────────────
namespace {

constexpr double kPi = 3.14159265358979323846;

// Logistic squash to [0,1] centred at x0 with steepness k.
inline double sigmoid(double x, double x0, double k) {
    return 1.0 / (1.0 + std::exp(-k * (x - x0)));
}

inline double clamp01(double v) {
    return std::max(0.0, std::min(1.0, v));
}

// Minimum-image displacement along one axis.
inline double minImage(double d, double L, bool periodic) {
    if (!periodic || L <= 0.0) return d;
    return d - L * std::round(d / L);
}

inline double pbcDist2(const std::array<double,3>& a,
                       const std::array<double,3>& b,
                       const SimBox& box) {
    double dx = minImage(a[0] - b[0], box.lx(), box.periodic[0]);
    double dy = minImage(a[1] - b[1], box.ly(), box.periodic[1]);
    double dz = minImage(a[2] - b[2], box.lz(), box.periodic[2]);
    return dx*dx + dy*dy + dz*dz;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// HybridParams::applyPreset
// ─────────────────────────────────────────────────────────────────────────────
void HybridParams::applyPreset(const std::string& name) {
    // All presets overwrite every weight; the cloud-context signals are only
    // enabled by 'survival' so the historical presets keep their invariants.
    w_centrality = 0.0; w_cloud_density = 0.0; ws_auto_accept = true;

    if (name == "ws") {
        w_ws = 1.0; w_soft_ws = 0.0; w_voronoi = 0.0; w_density = 0.0;
        w_soap = 0.0; w_topology = 0.0; w_transit = 0.0; w_frenkel = 0.0;
        accept_threshold = 0.5;
    } else if (name == "robust") {
        w_ws = 1.0; w_soft_ws = 0.8; w_voronoi = 0.6; w_density = 0.6;
        w_soap = 0.5; w_topology = 0.7; w_transit = 1.5; w_frenkel = 0.3;
        accept_threshold = 0.5;
    } else if (name == "relaxed") {
        // For cooled / relaxed final structures (e.g. post-cascade finalCool):
        // the MD has already resolved which Frenkel pairs recombined, so the
        // surviving vacancies are stable real defects. We keep the thermal-noise
        // robust positive signals but DISABLE the ballistic transit and Frenkel
        // recombination penalties (which are only meaningful for a ballistic /
        // peak-damage snapshot). This reconciles the hybrid count with WS while
        // staying more robust to thermal noise than raw WS.
        w_ws = 1.0; w_soft_ws = 0.8; w_voronoi = 0.6; w_density = 0.6;
        w_soap = 0.5; w_topology = 0.7; w_transit = 0.0; w_frenkel = 0.0;
        accept_threshold = 0.5;
    } else if (name == "sensitive") {
        w_ws = 1.0; w_soft_ws = 1.0; w_voronoi = 0.8; w_density = 1.0;
        w_soap = 0.7; w_topology = 0.8; w_transit = 0.8; w_frenkel = 0.2;
        accept_threshold = 0.35;
    } else if (name == "survival") {
        // Rank peak-damage candidates by survival likelihood. Weights follow
        // what per-defect survival labels supported (tracked FeCrNi 4–8 keV
        // cascades, leave-one-cascade-out AUC ≈ 0.70): the dominant predictor
        // is cloud centrality (core-shell structure), then local candidate
        // density and the KDE density deficit. The saturated signals
        // (soft_ws/voronoi/soap/topology) and both penalties carried no
        // usable peak-frame information and stay off. ws_auto_accept is
        // disabled so the reported count reflects the score threshold —
        // i.e. the PREDICTED SURVIVORS, not every WS-vacant site.
        w_ws = 0.3; w_soft_ws = 0.0; w_voronoi = 0.0; w_density = 0.4;
        w_soap = 0.0; w_topology = 0.0; w_transit = 0.0; w_frenkel = 0.0;
        w_centrality = 1.0; w_cloud_density = 0.4;
        accept_threshold = 0.6;
        ws_auto_accept = false;
    } else {
        throw std::invalid_argument(
            "HybridParams::applyPreset: unknown preset '" + name +
            "' (valid: ws, robust, relaxed, sensitive, survival)");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Acceptance rule
// ─────────────────────────────────────────────────────────────────────────────
//   A candidate is accepted iff:
//     (a) no transit / Frenkel penalty is firing, AND
//     (b) the consensus score passes the threshold, OR ws_auto_accept is set
//         and the site is WS-vacant (preserves the WS-equivalence invariants
//         of the ws/robust/relaxed presets; 'survival' disables it so the
//         count reflects the predicted survivors).
static inline bool isAccepted(const HybridVacancy& c, const HybridParams& p) {
    const bool penalty_clear = !c.transit_filtered && !c.in_recombination;
    const bool ws_says       = p.ws_auto_accept && c.breakdown.ws >= 0.5;
    const bool score_passes  = c.consensus_score >= p.accept_threshold;
    return penalty_clear && (ws_says || score_passes);
}

// ─────────────────────────────────────────────────────────────────────────────
// HybridVacancyDetector
// ─────────────────────────────────────────────────────────────────────────────
HybridVacancyDetector::HybridVacancyDetector(const HybridParams& p)
    : p_(p) {}

// ─── buildReference ──────────────────────────────────────────────────────────
void HybridVacancyDetector::buildReference(const Frame& ref_frame,
                                            const WignerSeitz& ws_built,
                                            const ReferenceSet* soap_ref) {
    if (ref_frame.size() == 0)
        throw std::invalid_argument("HybridVacancyDetector::buildReference: empty frame");

    ref_      = &ref_frame;
    ws_       = &ws_built;
    soap_ref_ = soap_ref;

    const SimBox& box = ref_frame.box;
    const double Lx = box.lx(), Ly = box.ly(), Lz = box.lz();
    const double V  = Lx * Ly * Lz;
    const int    N  = ref_frame.size();

    vol_per_atom_ = (N > 0) ? V / N : 0.0;

    // Heuristic 1st-neighbour distance: rough estimate from density,
    // then refined by the median nearest-neighbour distance of a sample.
    // For BCC: a ≈ (2V/N)^{1/3}, NN = a·√3/2.
    // For FCC: a ≈ (4V/N)^{1/3}, NN = a/√2.
    // We use BCC as default but the median sample is what actually drives nn_ref_.
    const double a_bcc = std::cbrt(2.0 * vol_per_atom_);
    const double nn_guess = a_bcc * std::sqrt(3.0) / 2.0;

    // Build a CellList over the reference frame at ~1.6·nn_guess to cover the
    // 1st shell with margin, then collect the minimum distance per atom.
    const double r_search = std::max(1.6 * nn_guess, 3.0);
    CellList cl;
    cl.build(ref_frame, r_search);

    // Sample up to 2000 atoms to keep this O(2000 · neighbours).
    const int sample_max = std::min(N, 2000);
    const int stride     = std::max(1, N / sample_max);

    std::vector<double> nn_dists;
    nn_dists.reserve(sample_max);
    std::vector<std::array<double,4>> tmp;

    for (int i = 0; i < N; i += stride) {
        const auto& a = ref_frame.atoms[i];
        tmp.clear();
        cl.neighborsInRadius(a.x, a.y, a.z, r_search, tmp);
        double min_r2 = std::numeric_limits<double>::max();
        for (const auto& v : tmp) {
            if (v[3] > 1e-10 && v[3] < min_r2) min_r2 = v[3];   // skip self (r²≈0)
        }
        if (min_r2 < std::numeric_limits<double>::max())
            nn_dists.push_back(std::sqrt(min_r2));
    }

    if (nn_dists.empty()) {
        nn_ref_ = nn_guess;
    } else {
        std::sort(nn_dists.begin(), nn_dists.end());
        nn_ref_ = nn_dists[nn_dists.size() / 2];     // median
    }

    // Mean coordination number within 1.2·nn_ref_ (BCC ⇒ 8, FCC ⇒ 12).
    double coord_sum = 0.0;
    int    coord_n   = 0;
    const double r_shell = 1.2 * nn_ref_;
    for (int i = 0; i < N; i += stride) {
        const auto& a = ref_frame.atoms[i];
        tmp.clear();
        cl.neighborsInRadius(a.x, a.y, a.z, r_shell, tmp);
        int cnt = 0;
        for (const auto& v : tmp) if (v[3] > 1e-10) ++cnt;       // skip self
        coord_sum += cnt;
        ++coord_n;
    }
    coord_ref_ = (coord_n > 0) ? coord_sum / coord_n : 8.0;

    // σ_thermal: explicit override > soap_ref-derived > 0.05·nn_ref default.
    if (p_.thermal_sigma > 0.0) {
        sigma_thermal_ = p_.thermal_sigma;
    } else {
        sigma_thermal_ = 0.05 * nn_ref_;
    }
    // Guard against pathological values
    if (sigma_thermal_ <= 1e-6) sigma_thermal_ = 0.05 * nn_ref_;
}

// ─── detect ──────────────────────────────────────────────────────────────────
std::vector<HybridVacancy> HybridVacancyDetector::detect(
    const Frame& dmg_frame,
    const std::vector<WSAtomResult>& ws_results)
{
    if (!ws_ || !ref_)
        throw std::runtime_error("HybridVacancyDetector::detect: buildReference() not called");

    cand_.clear();

    const SimBox& box = dmg_frame.box;
    const double sigma2 = 2.0 * sigma_thermal_ * sigma_thermal_;
    const double h      = std::max(0.25 * nn_ref_, 0.5);   // KDE bandwidth
    const double h2     = 2.0 * h * h;
    const double kde_norm = 1.0 / (std::pow(2.0 * kPi, 1.5) * h * h * h);
    const double rho_ref  = (vol_per_atom_ > 0.0) ? 1.0 / vol_per_atom_ : 1.0;

    const double r_signal = std::max({3.0 * sigma_thermal_, 2.0 * nn_ref_, 3.0 * h, 1.2 * nn_ref_});

    // Build a CellList over damaged frame, with cutoff covering all signal radii.
    CellList cl_dmg;
    cl_dmg.build(dmg_frame, r_signal);

    // ── Step 1: assemble candidate list ──────────────────────────────────────
    //   (a) each WSSite with occupancy == 0
    //   (b) plus grid-vacancy clusters that do NOT coincide with a vacant WS site
    const auto& sites = ws_->sites();

    std::vector<HybridVacancy> raw;
    raw.reserve(sites.size() / 16 + 64);

    for (size_t s = 0; s < sites.size(); ++s) {
        if (sites[s].occupancy == 0) {
            HybridVacancy hv;
            hv.pos[0] = sites[s].x;
            hv.pos[1] = sites[s].y;
            hv.pos[2] = sites[s].z;
            hv.ref_site_idx = static_cast<int>(s);
            raw.push_back(hv);
        }
    }

    // List of vacant WS sites (used both for candidate dedup and for the
    // transit-penalty: an atom near a vacant site is "in transit").
    const double dedup_r2 = (0.5 * nn_ref_) * (0.5 * nn_ref_);

    // Optionally seed with grid-vacancy clusters that miss any vacant WS site.
    //
    // The grid threshold must be > NN distance: otherwise EVERY octahedral /
    // tetrahedral interstitial hole in FCC and BCC lattices flags as a void
    // (~10⁶ false positives on a 1M-atom system).  We use 0.9·nn_ref, which
    // only fires on real missing-atom pockets.
    {
        DefectClassifier dc;
        auto grid_pts = dc.findVacanciesGrid(dmg_frame, p_.grid_spacing, 0.9 * nn_ref_);
        auto clusters = dc.clusterVacancyPoints(grid_pts, 0.9 * nn_ref_, box);

        // Spatial lookup of vacant WS sites so dedup is O(M_grid · 27) instead of O(M·N_ws).
        Frame ws_vac_frame;
        ws_vac_frame.box = box;
        for (const auto& hv : raw) {
            Atom a; a.id = 0; a.type = 0;
            a.x = hv.pos[0]; a.y = hv.pos[1]; a.z = hv.pos[2];
            ws_vac_frame.atoms.push_back(a);
        }
        const bool have_ws_vacs = !ws_vac_frame.atoms.empty();
        CellList cl_wsvac;
        if (have_ws_vacs) cl_wsvac.build(ws_vac_frame, std::max(nn_ref_, 1.0));

        for (const auto& cl : clusters) {
            if (have_ws_vacs) {
                const double d2 = cl_wsvac.nearestDist2FromPoint(
                    cl.center[0], cl.center[1], cl.center[2]);
                if (d2 <= dedup_r2) continue;
            }
            HybridVacancy hv;
            hv.pos = cl.center;
            hv.ref_site_idx = -1;
            raw.push_back(hv);
        }
    }

    // Precompute interstitial WS-site positions (occ ≥ 2) for Frenkel penalty.
    std::vector<std::array<double,3>> interstitials;
    interstitials.reserve(64);
    for (const auto& s : sites)
        if (s.occupancy >= 2) interstitials.push_back({s.x, s.y, s.z});

    // ── Step 2: per-candidate signals ────────────────────────────────────────
    std::vector<std::array<double,4>> nbrs;
    nbrs.reserve(128);

    // For the transit penalty we need to know whether an atom near a vacant
    // site is itself displaced (assigned to a different WS site that is now
    // crowded).  Pre-compute: is the WS site of atom j crowded (occ ≥ 2)?
    const bool have_ws_results = !ws_results.empty();

    for (auto& hv : raw) {
        const double qx = hv.pos[0], qy = hv.pos[1], qz = hv.pos[2];

        nbrs.clear();
        cl_dmg.neighborsInRadius(qx, qy, qz, r_signal, nbrs);

        // ── (1) WS binary
        hv.breakdown.ws = (hv.ref_site_idx >= 0) ? 1.0 : 0.0;

        // ── (2) Soft-WS gaussian occupancy
        //   Σ_j exp(−r²/(2σ_th²)) divided by the same sum for a "perfect" pristine
        //   environment (sum over a uniform crystal at distance r contributes
        //   approximately ρ_ref · ∫ exp(−r²/2σ²) 4πr² dr = ρ_ref·(2πσ²)^{3/2}).
        double soft_sum = 0.0;
        for (const auto& v : nbrs) soft_sum += std::exp(-v[3] / sigma2);
        const double soft_norm = rho_ref * std::pow(2.0 * kPi * sigma_thermal_ * sigma_thermal_, 1.5);
        const double omega = (soft_norm > 0.0) ? soft_sum / soft_norm : 0.0;
        hv.breakdown.soft_ws = clamp01(1.0 - omega);

        // ── (3) NN-excess (Voronoi proxy)
        //   For each atom within 2·nn_ref of the candidate, find its own nearest
        //   neighbour distance; take the median; compare against nn_ref.
        std::vector<double> nn_local;
        nn_local.reserve(16);
        const double r_local = 2.0 * nn_ref_;
        for (const auto& v : nbrs) {
            if (v[3] > r_local * r_local) continue;
            // atom at qx+v[0], qy+v[1], qz+v[2]
            std::vector<std::array<double,4>> nn2;
            cl_dmg.neighborsInRadius(qx + v[0], qy + v[1], qz + v[2], 1.8 * nn_ref_, nn2);
            double m2 = std::numeric_limits<double>::max();
            for (const auto& w : nn2)
                if (w[3] > 1e-10 && w[3] < m2) m2 = w[3];
            if (m2 < std::numeric_limits<double>::max())
                nn_local.push_back(std::sqrt(m2));
        }
        double vor_anom = 0.0;
        if (!nn_local.empty()) {
            std::sort(nn_local.begin(), nn_local.end());
            const double med = nn_local[nn_local.size() / 2];
            vor_anom = clamp01((med - nn_ref_) / nn_ref_);
        }
        hv.breakdown.voronoi_anomaly = vor_anom;

        // ── (4) Density deficit (KDE)
        double dens = 0.0;
        for (const auto& v : nbrs) {
            if (v[3] > (3.0 * h) * (3.0 * h)) continue;
            dens += std::exp(-v[3] / h2);
        }
        dens *= kde_norm;
        hv.breakdown.density_deficit = clamp01(1.0 - dens / rho_ref);

        // ── (5) SOAP neighbor anomaly
        //   Mean dist_to_ref over neighbours within 1.2·nn_ref. Converted to a
        //   rank percentil within the frame's candidate set in Step 2b. When
        //   DVs are not computed (dist_to_ref == 0 everywhere) every candidate
        //   ties and the signal collapses to the neutral 0.5.
        double soap_mean = 0.0;
        int    soap_n    = 0;
        std::vector<int>    nbr_idx;
        std::vector<double> nbr_r2;
        cl_dmg.neighborsInRadiusIndexed(qx, qy, qz, 1.2 * nn_ref_, nbr_idx, nbr_r2);
        for (int j : nbr_idx) {
            soap_mean += dmg_frame.atoms[j].dist_to_ref;
            ++soap_n;
        }
        if (soap_n > 0) soap_mean /= soap_n;

        // Store the RAW neighbourhood mean distance here; it is converted to a
        // rank percentil within this frame's candidate set in Step 2b. The old
        // absolute normalisation (reference mean + 3σ) saturated at 1.0 for
        // EVERY candidate in damage-peak frames (measured AUC 0.500 against
        // tracked survival labels in all 7 FeCrNi cascades) — the scale is
        // calibrated for near-lattice discrimination, not the cascade core.
        hv.breakdown.soap_neighbor = soap_mean;

        // ── (6) Topology / coordination anomaly
        //   For each atom within 1.5·nn_ref of the candidate, count its own neighbours
        //   within 1.2·nn_ref and compute |n - coord_ref| / coord_ref.
        double topo_sum = 0.0;
        int    topo_n   = 0;
        const double r_topo_query = 1.5 * nn_ref_;
        for (const auto& v : nbrs) {
            if (v[3] > r_topo_query * r_topo_query) continue;
            const double ax = qx + v[0], ay = qy + v[1], az = qz + v[2];
            std::vector<std::array<double,4>> nn2;
            cl_dmg.neighborsInRadius(ax, ay, az, 1.2 * nn_ref_, nn2);
            int cnt = 0;
            for (const auto& w : nn2) if (w[3] > 1e-10) ++cnt;
            topo_sum += std::abs(static_cast<double>(cnt) - coord_ref_) / coord_ref_;
            ++topo_n;
        }
        hv.breakdown.topology_anomaly = (topo_n > 0) ? clamp01(topo_sum / topo_n) : 0.0;

        // ── (P1) Transit penalty
        //   Fires only when there is an atom in the band [0.5·nn, transit_factor·nn]
        //   from the candidate AND that atom is itself displaced (either flagged
        //   as a WS interstitial — i.e. its WS cell has occ ≥ 2 — or has high SOAP
        //   defect_prob).  Pure lattice neighbours of a real vacancy do NOT trigger
        //   this penalty, so true vacancies in a cold lattice are preserved.
        double min_r2_band = std::numeric_limits<double>::max();
        const double r_lo = 0.5 * nn_ref_;
        const double r_hi = p_.transit_factor * nn_ref_;
        const double r_lo2 = r_lo * r_lo, r_hi2 = r_hi * r_hi;
        {
            std::vector<int>    idx_band;
            std::vector<double> r2_band;
            cl_dmg.neighborsInRadiusIndexed(qx, qy, qz, r_hi, idx_band, r2_band);
            for (size_t k = 0; k < idx_band.size(); ++k) {
                const double r2 = r2_band[k];
                if (r2 < r_lo2 || r2 > r_hi2) continue;
                const int j = idx_band[k];
                // A neighbour atom is "in transit" iff its WS cell is crowded
                // (occ ≥ 2, i.e. an interstitial).  defect_prob alone is NOT
                // sufficient: atoms first-shell-adjacent to a real vacancy
                // routinely have high defect_prob (their local environment is
                // distorted by the missing atom) but they are NOT in transit.
                bool displaced = false;
                if (have_ws_results) {
                    const auto& res = ws_results[j];
                    if (res.is_interstitial) displaced = true;
                }
                if (!displaced) continue;
                if (r2 < min_r2_band) min_r2_band = r2;
            }
        }
        if (min_r2_band < std::numeric_limits<double>::max()) {
            const double r = std::sqrt(min_r2_band);
            hv.breakdown.transit_penalty = clamp01((r_hi - r) / (r_hi - r_lo));
        } else {
            hv.breakdown.transit_penalty = 0.0;
        }

        // ── (P2) Frenkel coherence
        double min_inter_r2 = std::numeric_limits<double>::max();
        for (const auto& ip : interstitials) {
            const double r2 = pbcDist2(hv.pos, ip, box);
            if (r2 < min_inter_r2) min_inter_r2 = r2;
        }
        const double rrec2 = p_.recomb_radius * p_.recomb_radius;
        hv.breakdown.frenkel_penalty =
            (min_inter_r2 < std::numeric_limits<double>::max())
                ? std::exp(-min_inter_r2 / (2.0 * rrec2))
                : 0.0;

    }

    // ── Step 2b: cloud-context signals (cascade core-shell structure) ───────
    // Reference cloud = WS-vacant candidates (the population whose survival
    // was tracked). Signals are rank percentiles WITHIN this frame: the
    // operative question is relative ("which of THESE candidates survive"),
    // and absolute distances shift with PKA energy.
    {
        std::vector<int> ws_idx;
        ws_idx.reserve(raw.size());
        for (size_t i = 0; i < raw.size(); ++i)
            if (raw[i].ref_site_idx >= 0) ws_idx.push_back(static_cast<int>(i));

        const int n_all = static_cast<int>(raw.size());

        // Rank percentil con empates promediados (como pandas rank pct).
        // Columnas constantes (sin información) quedan en 0.5 neutro.
        auto rankPct = [n_all](const std::vector<double>& v) {
            std::vector<int> order(n_all);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](int a, int b){ return v[a] < v[b]; });
            std::vector<double> pct(n_all, 0.5);
            if (v[order.front()] == v[order.back()])
                return pct;                          // todo empatado: neutro
            int i = 0;
            while (i < n_all) {
                int j = i;
                while (j + 1 < n_all && v[order[j+1]] == v[order[i]]) ++j;
                const double avg_rank = 0.5 * (i + j) + 1.0;   // 1-based
                for (int k = i; k <= j; ++k)
                    pct[order[k]] = avg_rank / n_all;
                i = j + 1;
            }
            return pct;
        };

        // Señal 5 (SOAP): rank del promedio crudo de dist_to_ref vecinal.
        if (n_all > 1) {
            std::vector<double> soap_raw(n_all);
            for (int i = 0; i < n_all; ++i)
                soap_raw[i] = raw[i].breakdown.soap_neighbor;
            const auto pct_s = rankPct(soap_raw);
            for (int i = 0; i < n_all; ++i)
                raw[i].breakdown.soap_neighbor = pct_s[i];
        } else if (n_all == 1) {
            raw[0].breakdown.soap_neighbor = 0.5;
        }

        if (!ws_idx.empty() && n_all > 1) {
            // Centroid of the WS-vacancy cloud under PBC via the circular
            // mean per axis (robust regardless of where the cloud sits in
            // the box — a single-reference unwrap folds clouds that span
            // more than half a box length).
            const double L[3]  = { box.lx(), box.ly(), box.lz() };
            const double lo[3] = { box.xb[0], box.yb[0], box.zb[0] };
            std::array<double,3> cent{0,0,0};
            for (int d = 0; d < 3; ++d) {
                if (box.periodic[d] && L[d] > 0.0) {
                    double sc = 0.0, ss = 0.0;
                    for (int i : ws_idx) {
                        const double th = 2.0 * kPi * (raw[i].pos[d] - lo[d]) / L[d];
                        sc += std::cos(th);
                        ss += std::sin(th);
                    }
                    double th = std::atan2(ss, sc);
                    if (th < 0.0) th += 2.0 * kPi;
                    cent[d] = lo[d] + th * L[d] / (2.0 * kPi);
                } else {
                    double s = 0.0;
                    for (int i : ws_idx) s += raw[i].pos[d];
                    cent[d] = s / ws_idx.size();
                }
            }

            // Distances with plain minimum image — no unwrap needed.
            std::vector<double> dcent(n_all), ldens(n_all);
            const double cr2 = p_.cloud_radius * p_.cloud_radius;
            for (int i = 0; i < n_all; ++i) {
                dcent[i] = std::sqrt(pbcDist2(raw[i].pos, cent, box));
                int cnt = 0;
                for (int j : ws_idx) {
                    if (j == i) continue;
                    if (pbcDist2(raw[i].pos, raw[j].pos, box) <= cr2) ++cnt;
                }
                ldens[i] = cnt;
            }

            const auto pct_d = rankPct(dcent);
            const auto pct_l = rankPct(ldens);
            for (int i = 0; i < n_all; ++i) {
                raw[i].breakdown.cloud_centrality = 1.0 - pct_d[i];
                raw[i].breakdown.cloud_density    = pct_l[i];
            }
        }
    }

    // ── Step 3: consensus + accept decision ──────────────────────────────────
    const double w_sum =
          p_.w_ws + p_.w_soft_ws + p_.w_voronoi + p_.w_density
        + p_.w_soap + p_.w_topology
        + p_.w_centrality + p_.w_cloud_density;

    for (auto& hv : raw) {
        double pos_score = 0.0;
        if (w_sum > 0.0) {
            pos_score = (
                  p_.w_ws            * hv.breakdown.ws
                + p_.w_soft_ws       * hv.breakdown.soft_ws
                + p_.w_voronoi       * hv.breakdown.voronoi_anomaly
                + p_.w_density       * hv.breakdown.density_deficit
                + p_.w_soap          * hv.breakdown.soap_neighbor
                + p_.w_topology      * hv.breakdown.topology_anomaly
                + p_.w_centrality    * hv.breakdown.cloud_centrality
                + p_.w_cloud_density * hv.breakdown.cloud_density
            ) / w_sum;
        }
        const double penalty = p_.w_transit * hv.breakdown.transit_penalty
                             + p_.w_frenkel * hv.breakdown.frenkel_penalty;

        hv.consensus_score          = pos_score - penalty;
        hv.breakdown.consensus_score = hv.consensus_score;

        hv.transit_filtered  = (p_.w_transit > 0.0) && (hv.breakdown.transit_penalty > 0.5);
        hv.in_recombination  = (p_.w_frenkel > 0.0) && (hv.breakdown.frenkel_penalty > 0.5);
        hv.accepted          = isAccepted(hv, p_);
    }

    // ── Step 4: sort by score ────────────────────────────────────────────────
    std::sort(raw.begin(), raw.end(), [](const HybridVacancy& a, const HybridVacancy& b) {
        return a.consensus_score > b.consensus_score;
    });

    cand_ = std::move(raw);
    return cand_;
}

// ─────────────────────────────────────────────────────────────────────────────
// Accessors — counts use the accept decision stored by detect().
// ─────────────────────────────────────────────────────────────────────────────
int HybridVacancyDetector::vacancyCount() const {
    int n = 0;
    for (const auto& c : cand_) if (c.accepted) ++n;
    return n;
}

int HybridVacancyDetector::wsAgreeCount() const {
    int n = 0;
    for (const auto& c : cand_)
        if (c.ref_site_idx >= 0 && c.accepted) ++n;
    return n;
}

int HybridVacancyDetector::wsOnlyCount() const {
    int n = 0;
    for (const auto& c : cand_)
        if (c.ref_site_idx >= 0 && !c.accepted) ++n;
    return n;
}

int HybridVacancyDetector::hybridOnlyCount() const {
    int n = 0;
    for (const auto& c : cand_)
        if (c.ref_site_idx < 0 && c.accepted) ++n;
    return n;
}

} // namespace DistTool
