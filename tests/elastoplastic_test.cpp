#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/ELASTOPLASTIC.h"

TEST_CASE("ELASTOPLASTIC::ComputePolarDecomposition splits F into a rotation R and symmetric S", "[elastoplastic]") {
    ELASTOPLASTIC ConstitutiveModel;

    Matrix2d sample;
    sample << 2.0, 0.5,
              0.3, 1.5;
    REQUIRE(sample.determinant() > 0);  // physically valid deformation gradients have J > 0

    vector<Matrix2d> F = {sample};
    vector<Matrix2d> R(F.size());
    vector<Matrix2d> S(F.size());

    ConstitutiveModel.ComputePolarDecomposition(F, R, S);

    // F == R * S
    REQUIRE((R[0] * S[0]).isApprox(F[0]));

    // R is a proper rotation: orthogonal (R * R^T == I) with det(R) == 1,
    // not just any orthogonal matrix (which could be a reflection).
    REQUIRE((R[0] * R[0].transpose()).isApprox(Matrix2d::Identity()));
    REQUIRE(R[0].determinant() == Catch::Approx(1.0));

    // S is symmetric
    REQUIRE(S[0].isApprox(S[0].transpose()));
}

// mu(J_P) = mu0 * exp(hardeningCoeff * (1 - J_P)), from Stomakhin et al. 2013
// (ACM TOG 32,4, Article 102), eq. 2.
TEST_CASE("ELASTOPLASTIC::ComputeMu applies exponential hardening based on J_P", "[elastoplastic]") {
    ELASTOPLASTIC ConstitutiveModel;
    const double mu0 = 100.0;
    const double hardeningCoeff = 10.0;
    vector<double> JP = {1.0, 0.95, 1.05};

    SECTION("plasticity disabled: mu is just mu0, JP is ignored") {
        vector<double> mu = ConstitutiveModel.ComputeMu(false, mu0, hardeningCoeff, JP);

        REQUIRE(mu.size() == JP.size());
        for (double value : mu) {
            REQUIRE(value == Catch::Approx(mu0));
        }
    }

    SECTION("plasticity enabled: mu follows the exponential hardening formula") {
        vector<double> mu = ConstitutiveModel.ComputeMu(true, mu0, hardeningCoeff, JP);

        REQUIRE(mu.size() == JP.size());
        REQUIRE(mu[0] == Catch::Approx(mu0));  // J_P == 1: no plastic deformation yet, no hardening
        REQUIRE(mu[1] == Catch::Approx(mu0 * exp(hardeningCoeff * (1.0 - JP[1]))));  // compressed: stiffer
        REQUIRE(mu[2] == Catch::Approx(mu0 * exp(hardeningCoeff * (1.0 - JP[2]))));  // stretched: softer
        REQUIRE(mu[1] > mu0);
        REQUIRE(mu[2] < mu0);
    }
}
