#include "WignerSeitz.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DistTool {

// ─────────────────────────────────────────────────────────────────────────────

void WignerSeitz::build(const Frame& ref_frame, double r_cut) {
    if (r_cut <= 0.0)
        throw std::runtime_error("WignerSeitz::build: r_cut must be > 0");

    r_cut_  = r_cut;
    ref_box_ = ref_frame.box;

    Lx_  = ref_box_.lx();  Ly_  = ref_box_.ly();  Lz_  = ref_box_.lz();
    xlo_ = ref_box_.xb[0]; ylo_ = ref_box_.yb[0]; zlo_ = ref_box_.zb[0];

    nx_ = std::max(1, static_cast<int>(std::floor(Lx_ / r_cut_)));
    ny_ = std::max(1, static_cast<int>(std::floor(Ly_ / r_cut_)));
    nz_ = std::max(1, static_cast<int>(std::floor(Lz_ / r_cut_)));
    cx_ = Lx_ / nx_;
    cy_ = Ly_ / ny_;
    cz_ = Lz_ / nz_;

    const int N = ref_frame.size();
    sites_.clear();
    sites_.reserve(N);
    cells_.assign(nx_ * ny_ * nz_, {});

    for (int i = 0; i < N; ++i) {
        const auto& a = ref_frame.atoms[i];
        sites_.push_back({a.id, a.x, a.y, a.z, 0, 0.0});

        int ix, iy, iz;
        posCell(a.x, a.y, a.z, ix, iy, iz);
        cells_[cellIdx(ix, iy, iz)].push_back(i);
    }
}

// ─────────────────────────────────────────────────────────────────────────────

std::vector<WSAtomResult> WignerSeitz::classify(const Frame& dmg_frame) {
    // Reset occupancies from a previous classify() call.
    for (auto& s : sites_) { s.occupancy = 0; s.min_dist = 0.0; }

    const int N = dmg_frame.size();
    std::vector<WSAtomResult> results(N);

    // Step 1: assign each damaged atom to its nearest reference site.
    for (int i = 0; i < N; ++i) {
        const auto& a = dmg_frame.atoms[i];
        double dist;
        int si = nearestSite(a.x, a.y, a.z, dist);

        results[i].ref_site_idx = si;
        results[i].dist         = dist;

        if (si >= 0) {
            ++sites_[si].occupancy;
            // Track distance of nearest occupant for reporting.
            if (sites_[si].occupancy == 1 || dist < sites_[si].min_dist)
                sites_[si].min_dist = dist;
        }
    }

    // Step 2: annotate each result with the final occupancy of its WS cell.
    for (int i = 0; i < N; ++i) {
        int si = results[i].ref_site_idx;
        if (si < 0) continue;
        results[i].ws_occ        = sites_[si].occupancy;
        results[i].is_interstitial = (sites_[si].occupancy >= 2);
    }

    return results;
}

// ─────────────────────────────────────────────────────────────────────────────

int WignerSeitz::vacancyCount() const {
    int n = 0;
    for (const auto& s : sites_)
        if (s.occupancy == 0) ++n;
    return n;
}

int WignerSeitz::interstitialCount() const {
    int n = 0;
    for (const auto& s : sites_)
        if (s.occupancy >= 2) n += s.occupancy - 1;
    return n;
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

int WignerSeitz::nearestSite(double qx, double qy, double qz,
                              double& dist_out) const {
    int ix, iy, iz;
    posCell(qx, qy, qz, ix, iy, iz);

    double min_r2  = std::numeric_limits<double>::max();
    int    nearest = -1;

    for (int dix = -1; dix <= 1; ++dix)
    for (int diy = -1; diy <= 1; ++diy)
    for (int diz = -1; diz <= 1; ++diz) {
        int jx = ix + dix, jy = iy + diy, jz = iz + diz;
        if (ref_box_.periodic[0]) jx = ((jx % nx_) + nx_) % nx_;
        else if (jx < 0 || jx >= nx_) continue;
        if (ref_box_.periodic[1]) jy = ((jy % ny_) + ny_) % ny_;
        else if (jy < 0 || jy >= ny_) continue;
        if (ref_box_.periodic[2]) jz = ((jz % nz_) + nz_) % nz_;
        else if (jz < 0 || jz >= nz_) continue;

        for (int k : cells_[cellIdx(jx, jy, jz)]) {
            double dx = sites_[k].x - qx;
            double dy = sites_[k].y - qy;
            double dz = sites_[k].z - qz;
            if (ref_box_.periodic[0] && Lx_ > 0.0) dx -= Lx_ * std::round(dx / Lx_);
            if (ref_box_.periodic[1] && Ly_ > 0.0) dy -= Ly_ * std::round(dy / Ly_);
            if (ref_box_.periodic[2] && Lz_ > 0.0) dz -= Lz_ * std::round(dz / Lz_);
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < min_r2) { min_r2 = r2; nearest = k; }
        }
    }

    // Fallback: linear scan when either no site was found in the 27-cell shell
    // OR the best candidate is farther than r_cut_ (the shell guarantee only holds
    // for sites within r_cut_; beyond that, the true NN might be in a farther cell).
    if (nearest < 0 || min_r2 > r_cut_ * r_cut_) {
        for (int k = 0; k < static_cast<int>(sites_.size()); ++k) {
            double dx = sites_[k].x - qx;
            double dy = sites_[k].y - qy;
            double dz = sites_[k].z - qz;
            if (ref_box_.periodic[0] && Lx_ > 0.0) dx -= Lx_ * std::round(dx / Lx_);
            if (ref_box_.periodic[1] && Ly_ > 0.0) dy -= Ly_ * std::round(dy / Ly_);
            if (ref_box_.periodic[2] && Lz_ > 0.0) dz -= Lz_ * std::round(dz / Lz_);
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < min_r2) { min_r2 = r2; nearest = k; }
        }
    }

    dist_out = (nearest >= 0) ? std::sqrt(min_r2) : std::numeric_limits<double>::max();
    return nearest;
}

void WignerSeitz::posCell(double x, double y, double z,
                           int& ix, int& iy, int& iz) const {
    ix = static_cast<int>((x - xlo_) / cx_);
    iy = static_cast<int>((y - ylo_) / cy_);
    iz = static_cast<int>((z - zlo_) / cz_);
    ix = std::max(0, std::min(nx_ - 1, ix));
    iy = std::max(0, std::min(ny_ - 1, iy));
    iz = std::max(0, std::min(nz_ - 1, iz));
}

} // namespace DistTool
