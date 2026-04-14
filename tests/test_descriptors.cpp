#include "test_harness.h"
#include "RadialBasis.h"
#include "SOAPDescriptor.h"
#include "SphericalHarmonics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using namespace DistTool;

// ── SphericalHarmonics ────────────────────────────────────────────────────────

TEST_CASE("SphericalHarmonics::Y_0^0 is constant for all directions") {
    // Y_0^0 = 1/sqrt(4*pi) regardless of direction
    SphericalHarmonics sph(4);
    const double Y00 = 1.0 / std::sqrt(4.0 * M_PI);
    std::vector<double> out(sph.nComponents());

    for (auto& d : std::vector<std::array<double,3>>{
            {1,0,0}, {0,1,0}, {0,0,1},
            {1,1,0}, {1,1,1}, {2,-3,1}, {0.1, -0.5, 0.8}}) {
        sph.computeInto(d[0], d[1], d[2], out.data());
        APPROX_EQ(out[SphericalHarmonics::idx(0, 0)], Y00, 1e-12);
    }
}

TEST_CASE("SphericalHarmonics::Y_1 along z-axis") {
    // Along z-axis (sin(theta)=0): only m=0 components survive
    SphericalHarmonics sph(3);
    std::vector<double> out(sph.nComponents());
    sph.computeInto(0.0, 0.0, 1.0, out.data());

    // Y_1^0(z-hat) = sqrt(3/(4*pi))
    APPROX_EQ(out[SphericalHarmonics::idx(1,  0)],
              std::sqrt(3.0 / (4.0 * M_PI)), 1e-12);
    // Y_1^{+1} and Y_1^{-1} vanish (sin theta = 0)
    APPROX_EQ(out[SphericalHarmonics::idx(1,  1)], 0.0, 1e-12);
    APPROX_EQ(out[SphericalHarmonics::idx(1, -1)], 0.0, 1e-12);
}

TEST_CASE("SphericalHarmonics::addition theorem — sum_m Y_l^m^2 = (2l+1)/(4pi)") {
    // For real tesseral harmonics: sum_{m=-l}^{l} [Y_l^m(r)]^2 = (2l+1)/(4*pi)
    // for any unit direction r.
    SphericalHarmonics sph(6);
    std::vector<double> out(sph.nComponents());

    for (auto& d : std::vector<std::array<double,3>>{
            {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}, {2,-1,3}}) {
        sph.computeInto(d[0], d[1], d[2], out.data());
        for (int l = 0; l <= 6; ++l) {
            double sum = 0.0;
            for (int m = -l; m <= l; ++m)
                sum += out[SphericalHarmonics::idx(l, m)] *
                       out[SphericalHarmonics::idx(l, m)];
            APPROX_EQ(sum, (2.0*l + 1.0) / (4.0*M_PI), 1e-11);
        }
    }
}

// ── RadialBasis ───────────────────────────────────────────────────────────────

TEST_CASE("RadialBasis::cutoff is 1 at r=0 and 0 at r=r_cut") {
    RadialBasis rb(5, 5.0);
    APPROX_EQ(rb.cutoff(0.0),  1.0, 1e-12);
    APPROX_EQ(rb.cutoff(5.0),  0.0, 1e-12);
    APPROX_EQ(rb.cutoff(10.0), 0.0, 1e-12);
}

TEST_CASE("RadialBasis::computeInto returns zeros and false for r >= r_cut") {
    RadialBasis rb(5, 5.0);
    std::vector<double> phi(5, 99.0);
    bool ok = rb.computeInto(5.0, phi.data());
    REQUIRE(!ok);
    for (double v : phi) APPROX_EQ(v, 0.0, 1e-12);
}

TEST_CASE("RadialBasis::centres are evenly spaced in (0, r_cut]") {
    const int n_max = 6;
    const double r_cut = 5.0;
    RadialBasis rb(n_max, r_cut);
    // Centre n = (n+1) * r_cut / (n_max+1)
    std::vector<double> phi(n_max);
    // At each centre r_n, the n-th basis function should be the largest
    for (int n = 0; n < n_max; ++n) {
        double r_n = (n + 1.0) * r_cut / (n_max + 1.0);
        rb.computeInto(r_n, phi.data());
        double max_val = *std::max_element(phi.begin(), phi.end());
        APPROX_EQ(phi[n], max_val, 1e-10);
    }
}

// ── SOAPDescriptor ────────────────────────────────────────────────────────────

// Returns all 8 BCC nearest-neighbour displacement vectors for Fe (a=2.85 Å).
static std::vector<std::array<double,3>> bcc_neighbours() {
    const double a = 1.425;  // a/2 for Fe BCC
    std::vector<std::array<double,3>> nb;
    for (int sx : {-1,1}) for (int sy : {-1,1}) for (int sz : {-1,1})
        nb.push_back({sx*a, sy*a, sz*a});
    return nb;
}

TEST_CASE("SOAPDescriptor::DV is unit-length after normalisation") {
    SOAPParams p;  p.n_max=5; p.l_max=4; p.r_cut=5.0;
    SOAPDescriptor desc(p);
    auto dv = desc.compute(bcc_neighbours());
    double norm2 = 0.0;
    for (double v : dv) norm2 += v*v;
    APPROX_EQ(norm2, 1.0, 1e-12);
}

TEST_CASE("SOAPDescriptor::empty neighbourhood gives zero DV") {
    SOAPParams p;  p.n_max=5; p.l_max=4; p.r_cut=5.0;
    SOAPDescriptor desc(p);
    auto dv = desc.compute({});
    for (double v : dv) APPROX_EQ(v, 0.0, 1e-12);
}

TEST_CASE("SOAPDescriptor::rotational invariance — 90-degree rotation around z") {
    // BCC NN set has 4-fold symmetry around z: rotating by 90° permutes the
    // vectors but does not change the set, so the SOAP descriptor must be identical.
    SOAPParams p;  p.n_max=5; p.l_max=5; p.r_cut=5.0;
    SOAPDescriptor desc(p);

    auto orig = bcc_neighbours();

    // Rotate each vector 90° around z: (x,y,z) -> (-y, x, z)
    auto rotated = orig;
    for (auto& v : rotated) {
        double x = v[0], y = v[1];
        v[0] = -y;
        v[1] =  x;
    }

    auto dv_orig    = desc.compute(orig);
    auto dv_rotated = desc.compute(rotated);

    REQUIRE(dv_orig.size() == dv_rotated.size());
    for (size_t i = 0; i < dv_orig.size(); ++i)
        APPROX_EQ(dv_orig[i], dv_rotated[i], 1e-10);
}

TEST_CASE("SOAPDescriptor::same environment gives distance zero") {
    SOAPParams p;  p.n_max=5; p.l_max=4; p.r_cut=5.0;
    SOAPDescriptor desc(p);
    auto nb  = bcc_neighbours();
    auto dv1 = desc.compute(nb);
    auto dv2 = desc.compute(nb);
    double dist = 0.0;
    for (size_t i = 0; i < dv1.size(); ++i) { double d = dv1[i]-dv2[i]; dist += d*d; }
    APPROX_EQ(dist, 0.0, 1e-20);
}
