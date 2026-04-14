#include "test_harness.h"
#include "Statistics.h"
#include <cmath>
#include <vector>

using namespace DistTool;

// ── mean ─────────────────────────────────────────────────────────────────────

TEST_CASE("Statistics::mean — basic average") {
    auto m = Statistics::mean({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE(m.size() == 2);
    APPROX_EQ(m[0], 2.0, 1e-12);
    APPROX_EQ(m[1], 3.0, 1e-12);
}

TEST_CASE("Statistics::mean — single vector returns itself") {
    auto m = Statistics::mean({{5.0, -1.0, 3.0}});
    APPROX_EQ(m[0],  5.0, 1e-12);
    APPROX_EQ(m[1], -1.0, 1e-12);
    APPROX_EQ(m[2],  3.0, 1e-12);
}

TEST_CASE("Statistics::mean — throws on inconsistent vector sizes") {
    CHECK_THROWS(Statistics::mean({{1.0, 2.0}, {3.0, 4.0, 5.0}}));
}

// ── euclidean ─────────────────────────────────────────────────────────────────

TEST_CASE("Statistics::euclidean — 3-4-5 triangle") {
    APPROX_EQ(Statistics::euclidean({3.0, 0.0}, {0.0, -4.0}), 5.0, 1e-12);
}

TEST_CASE("Statistics::euclidean — identical vectors give 0") {
    std::vector<double> v(50, 1.7);
    APPROX_EQ(Statistics::euclidean(v, v), 0.0, 1e-12);
}

// ── chiProbability ────────────────────────────────────────────────────────────

TEST_CASE("Statistics::chiProbability — equals 1 at mode (k > 2)") {
    // Mode at d_peak = sigma * sqrt(k - 2)
    for (double k : {3.0, 5.0, 10.0, 50.0})
        for (double sigma : {0.05, 0.1, 0.5}) {
            double d_peak = sigma * std::sqrt(k - 2.0);
            APPROX_EQ(Statistics::chiProbability(d_peak, k, sigma), 1.0, 1e-9);
        }
}

TEST_CASE("Statistics::chiProbability — monotone decreasing for k <= 2") {
    // k=2: P = exp(-d^2/(2*sigma^2)), strictly decreasing
    double sigma = 0.1;
    double prev = Statistics::chiProbability(0.0, 2.0, sigma);
    for (double d : {0.05, 0.1, 0.2, 0.5}) {
        double p = Statistics::chiProbability(d, 2.0, sigma);
        REQUIRE(p < prev);
        prev = p;
    }
}

TEST_CASE("Statistics::chiProbability — farther from mode gives lower probability") {
    double k = 8.0, sigma = 0.1;
    double d_peak = sigma * std::sqrt(k - 2.0);
    double p_mode  = Statistics::chiProbability(d_peak,         k, sigma);
    double p_close = Statistics::chiProbability(d_peak + 0.02,  k, sigma);
    double p_far   = Statistics::chiProbability(d_peak + 0.2,   k, sigma);
    REQUIRE(p_mode > p_close);
    REQUIRE(p_close > p_far);
}

// ── fitChiParams ──────────────────────────────────────────────────────────────

TEST_CASE("Statistics::fitChiParams — recovers k=8, sigma=0.1 from exact moments") {
    // Construct distances with exact moments for chi(k=8, sigma=0.1):
    //   E[d^2] = k*sigma^2 = 0.08
    //   Var(d^2) = 2*k*sigma^4 = 0.0016
    // Two values d1=sqrt(0.12) and d2=0.2 in equal proportion satisfy these exactly.
    const double d1 = std::sqrt(0.12);  // sqrt(0.12)
    const double d2 = 0.2;
    std::vector<double> dists(100);
    for (int i = 0;  i < 50;  ++i) dists[i]    = d1;
    for (int i = 50; i < 100; ++i) dists[i]    = d2;

    double k_fit, sigma_fit;
    Statistics::fitChiParams(dists, k_fit, sigma_fit);
    APPROX_EQ(k_fit,     8.0, 1e-9);
    APPROX_EQ(sigma_fit, 0.1, 1e-9);
}
