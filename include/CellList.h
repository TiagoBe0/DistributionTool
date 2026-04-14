#pragma once
#include "AtomData.h"
#include <cmath>
#include <limits>
#include <vector>

namespace DistTool {

/**
 * Spatial cell list (linked-cell algorithm) for O(N) neighbour lookup.
 *
 * The simulation box is divided into cubic cells of side ≥ r_cut.
 * For each atom the search visits only the 27 surrounding cells.
 *
 * Complexity:
 *   build()     O(N)
 *   neighbours  O(27 · ρ · r_cut³)  per atom  (≈ constant for uniform density)
 *
 * Periodic boundary conditions are handled with minimum-image convention.
 */
class CellList {
public:
    CellList() = default;

    // Build the cell list from all atoms in `frame` for cutoff `r_cut`.
    void build(const Frame& frame, double r_cut) {
        frame_   = &frame;
        r_cut_   = r_cut;
        r_cut2_  = r_cut * r_cut;

        const SimBox& box = frame.box;
        Lx_ = box.lx();  Ly_ = box.ly();  Lz_ = box.lz();
        xlo_ = box.xb[0]; ylo_ = box.yb[0]; zlo_ = box.zb[0];

        // Number of cells per dimension (at least 1)
        nx_ = std::max(1, static_cast<int>(std::floor(Lx_ / r_cut)));
        ny_ = std::max(1, static_cast<int>(std::floor(Ly_ / r_cut)));
        nz_ = std::max(1, static_cast<int>(std::floor(Lz_ / r_cut)));

        cx_ = Lx_ / nx_;
        cy_ = Ly_ / ny_;
        cz_ = Lz_ / nz_;

        const int ncells = nx_ * ny_ * nz_;
        cells_.assign(ncells, {});

        // Assign each atom to its cell
        for (int i = 0; i < frame.size(); ++i) {
            int ix, iy, iz;
            atomCell(frame.atoms[i], ix, iy, iz);
            cells_[cellIdx(ix, iy, iz)].push_back(i);
        }
    }

    // Return squared distance to nearest atom in this cell list from an arbitrary
    // query point (qx, qy, qz).  Uses minimum-image PBC.  Build with
    // r_cut >= search_radius so the 3×3×3 cell shell covers the full search volume.
    // Returns std::numeric_limits<double>::max() if the frame is empty.
    double nearestDist2FromPoint(double qx, double qy, double qz) const {
        int ix, iy, iz;
        posCell(qx, qy, qz, ix, iy, iz);

        const SimBox& box = frame_->box;
        double min_r2 = std::numeric_limits<double>::max();

        for (int dix = -1; dix <= 1; ++dix)
        for (int diy = -1; diy <= 1; ++diy)
        for (int diz = -1; diz <= 1; ++diz) {
            int jx = ix + dix, jy = iy + diy, jz = iz + diz;
            if (box.periodic[0]) jx = ((jx % nx_) + nx_) % nx_;
            else if (jx < 0 || jx >= nx_) continue;
            if (box.periodic[1]) jy = ((jy % ny_) + ny_) % ny_;
            else if (jy < 0 || jy >= ny_) continue;
            if (box.periodic[2]) jz = ((jz % nz_) + nz_) % nz_;
            else if (jz < 0 || jz >= nz_) continue;

            for (int j : cells_[cellIdx(jx, jy, jz)]) {
                const Atom& aj = frame_->atoms[j];
                double dx = aj.x - qx;
                double dy = aj.y - qy;
                double dz = aj.z - qz;
                if (box.periodic[0] && Lx_ > 0.0) dx -= Lx_ * std::round(dx / Lx_);
                if (box.periodic[1] && Ly_ > 0.0) dy -= Ly_ * std::round(dy / Ly_);
                if (box.periodic[2] && Lz_ > 0.0) dz -= Lz_ * std::round(dz / Lz_);
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < min_r2) min_r2 = r2;
            }
        }
        return min_r2;
    }

    // Return displacement vectors (dx, dy, dz) from atom[i] to all neighbours
    // within r_cut (with PBC minimum-image).
    std::vector<std::array<double,3>> neighbours(int i) const {
        std::vector<std::array<double,3>> result;
        result.reserve(64);

        const Atom& ai = frame_->atoms[i];
        const SimBox& box = frame_->box;

        int ix, iy, iz;
        atomCell(ai, ix, iy, iz);

        for (int dix = -1; dix <= 1; ++dix)
        for (int diy = -1; diy <= 1; ++diy)
        for (int diz = -1; diz <= 1; ++diz) {
            int jx = ix + dix, jy = iy + diy, jz = iz + diz;

            // Wrap periodic dimensions; skip out-of-range cells in non-periodic ones.
            if (box.periodic[0]) jx = ((jx % nx_) + nx_) % nx_;
            else if (jx < 0 || jx >= nx_) continue;
            if (box.periodic[1]) jy = ((jy % ny_) + ny_) % ny_;
            else if (jy < 0 || jy >= ny_) continue;
            if (box.periodic[2]) jz = ((jz % nz_) + nz_) % nz_;
            else if (jz < 0 || jz >= nz_) continue;

            for (int j : cells_[cellIdx(jx, jy, jz)]) {
                if (j == i) continue;
                const Atom& aj = frame_->atoms[j];

                double dx = aj.x - ai.x;
                double dy = aj.y - ai.y;
                double dz = aj.z - ai.z;

                // Minimum-image convention
                if (box.periodic[0]) dx -= Lx_ * std::round(dx / Lx_);
                if (box.periodic[1]) dy -= Ly_ * std::round(dy / Ly_);
                if (box.periodic[2]) dz -= Lz_ * std::round(dz / Lz_);

                if (dx*dx + dy*dy + dz*dz < r_cut2_)
                    result.push_back({dx, dy, dz});
            }
        }
        return result;
    }

private:
    const Frame* frame_ = nullptr;
    double r_cut_ = 0, r_cut2_ = 0;
    double Lx_, Ly_, Lz_;
    double xlo_, ylo_, zlo_;
    int    nx_, ny_, nz_;
    double cx_, cy_, cz_;

    std::vector<std::vector<int>> cells_;

    int cellIdx(int ix, int iy, int iz) const {
        return ix * ny_ * nz_ + iy * nz_ + iz;
    }

    void posCell(double x, double y, double z, int& ix, int& iy, int& iz) const {
        ix = static_cast<int>((x - xlo_) / cx_);
        iy = static_cast<int>((y - ylo_) / cy_);
        iz = static_cast<int>((z - zlo_) / cz_);
        ix = std::max(0, std::min(nx_-1, ix));
        iy = std::max(0, std::min(ny_-1, iy));
        iz = std::max(0, std::min(nz_-1, iz));
    }

    void atomCell(const Atom& a, int& ix, int& iy, int& iz) const {
        posCell(a.x, a.y, a.z, ix, iy, iz);
    }
};

} // namespace DistTool
