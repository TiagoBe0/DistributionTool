// Microbenchmark: incremental (point-by-point) vs range insertion into a
// periodic 3D Delaunay triangulation, for growing N.
//
// Why this exists: the reference-free Delaunay void detector (--delaunay-voids)
// once took >98 min and never finished on a 1.31M-atom cascade frame. The cause
// was NOT a logic bug but the interaction of two things — point-by-point
// insertion does a "walk" point location, and a LAMMPS dump arrives in
// spatially-structured (crystal raster) order, which degrades the walk to O(N)
// per point => O(N^2) total. CGAL's range insert spatially sorts the batch
// (Hilbert curve) first, so it is immune to input order.
//
// This benchmark reproduces that: with random points the gap is ~1.5x; with
// crystal-raster order it blows up (~190x at N=20k, and grows as O(N^2)).
//
// Build:  g++ -O3 -march=native -std=c++17 scripts/bench_delaunay_insert.cpp \
//             -o /tmp/bench_delaunay_insert -lgmp -lmpfr
// Run:    /tmp/bench_delaunay_insert [N1 N2 ...]      # default 20k..160k
//         BENCH_MODE=0 /tmp/bench_delaunay_insert      # uniform random order
//         BENCH_MODE=1 /tmp/bench_delaunay_insert      # crystal raster (default)
//
// NOTE: in crystal-raster mode the incremental column scales O(N^2); keep N
// small (<= ~40k) unless you want to wait minutes.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using GT  = CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>;
using PDT = CGAL::Periodic_3_Delaunay_triangulation_3<GT>;
using P   = GT::Point_3;

// mode 0 = uniform random; mode 1 = crystal lattice in raster order (z,y,x)
// with small thermal jitter — mimics a LAMMPS dump ordered by atom id.
static std::vector<P> makePoints(int n, double L, int mode) {
    std::mt19937 rng(12345);
    std::vector<P> v; v.reserve(n);
    if (mode == 0) {
        std::uniform_real_distribution<double> u(0.0, L);
        for (int i = 0; i < n; ++i) v.emplace_back(u(rng), u(rng), u(rng));
    } else {
        const int side = (int)std::ceil(std::cbrt((double)n));
        const double a = L / side;                 // lattice spacing
        std::normal_distribution<double> jit(0.0, 0.02 * a);  // thermal noise
        for (int iz = 0; iz < side && (int)v.size() < n; ++iz)
        for (int iy = 0; iy < side && (int)v.size() < n; ++iy)
        for (int ix = 0; ix < side && (int)v.size() < n; ++ix) {
            double x = std::fmod((ix + 0.5) * a + jit(rng) + L, L);
            double y = std::fmod((iy + 0.5) * a + jit(rng) + L, L);
            double z = std::fmod((iz + 0.5) * a + jit(rng) + L, L);
            v.emplace_back(x, y, z);
        }
    }
    return v;
}

template <class F>
static double timed(F&& f) {
    auto t0 = std::chrono::steady_clock::now();
    f();
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

int main(int argc, char** argv) {
    const double L = 50.0;
    std::vector<int> Ns = {20000, 40000, 80000, 160000};
    if (argc > 1) { Ns.clear(); for (int i=1;i<argc;++i) Ns.push_back(std::atoi(argv[i])); }

    const int mode = (std::getenv("BENCH_MODE") ? std::atoi(std::getenv("BENCH_MODE")) : 1);
    std::printf("# point order: %s\n", mode == 0 ? "uniform random" : "crystal raster (z,y,x)");
    std::printf("%8s | %14s %7s | %14s %7s | %6s\n",
                "N", "incremental(ms)", "1cov", "range(ms)", "1cov", "ratio");
    for (int N : Ns) {
        auto pts = makePoints(N, L, mode);
        GT::Iso_cuboid_3 dom(0,0,0, L,L,L);

        bool inc_1cov = false, rng_1cov = false;
        double t_inc = timed([&]{
            PDT dt(dom);
            for (const auto& p : pts) dt.insert(p);
            inc_1cov = dt.is_1_cover();
        });
        double t_rng = timed([&]{
            PDT dt(dom);
            dt.insert(pts.begin(), pts.end(), true);
            rng_1cov = dt.is_1_cover();
        });
        std::printf("%8d | %14.1f %7d | %14.1f %7d | %5.1fx\n",
                    N, t_inc, (int)inc_1cov, t_rng, (int)rng_1cov, t_inc / t_rng);
        std::fflush(stdout);
    }
    return 0;
}
