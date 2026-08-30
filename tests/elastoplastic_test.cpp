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

// lambda(J_P) = lambda0 * exp(hardeningCoeff * (1 - J_P)), from Stomakhin
// et al. 2013 (ACM TOG 32,4, Article 102), eq. 2 — same shape as mu(J_P),
// applied to the other Lame parameter.
TEST_CASE("ELASTOPLASTIC::ComputeLambda applies exponential hardening based on J_P", "[elastoplastic]") {
    ELASTOPLASTIC ConstitutiveModel;
    const double lambda0 = 200.0;
    const double hardeningCoeff = 10.0;
    vector<double> JP = {1.0, 0.95, 1.05};

    SECTION("plasticity disabled: lambda is just lambda0, JP is ignored") {
        vector<double> lambda = ConstitutiveModel.ComputeLambda(false, lambda0, hardeningCoeff, JP);

        REQUIRE(lambda.size() == JP.size());
        for (double value : lambda) {
            REQUIRE(value == Catch::Approx(lambda0));
        }
    }

    SECTION("plasticity enabled: lambda follows the exponential hardening formula") {
        vector<double> lambda = ConstitutiveModel.ComputeLambda(true, lambda0, hardeningCoeff, JP);

        REQUIRE(lambda.size() == JP.size());
        REQUIRE(lambda[0] == Catch::Approx(lambda0));  // J_P == 1: no plastic deformation yet, no hardening
        REQUIRE(lambda[1] == Catch::Approx(lambda0 * exp(hardeningCoeff * (1.0 - JP[1]))));  // compressed: stiffer
        REQUIRE(lambda[2] == Catch::Approx(lambda0 * exp(hardeningCoeff * (1.0 - JP[2]))));  // stretched: softer
        REQUIRE(lambda[1] > lambda0);
        REQUIRE(lambda[2] < lambda0);
    }
}

// sigma = (2*mu/J) * (F_E - R) * F_E^T + (lambda/J) * (J_E - 1) * J_E * I,
// the fixed-corotated Cauchy stress from Stomakhin et al. 2013, eq. 6. At
// F_E = I, the polar decomposition gives R = I (so F_E - R = 0) and
// J_E = det(F_E) = 1 (so J_E - 1 = 0), so both terms vanish regardless of
// mu/lambda.
TEST_CASE("ELASTOPLASTIC::CauchyStress is zero at the identity deformation gradient", "[elastoplastic]") {
    ELASTOPLASTIC ConstitutiveModel;
    const double mu0 = 100.0;
    const double lambda0 = 200.0;
    const double hardeningCoeff = 10.0;

    vector<double> JElastic = {1.0};
    vector<double> JPlastic = {1.0};
    vector<Matrix2d> elasticDeformationGradient = {Matrix2d::Identity()};
    vector<Matrix2d> R(1);
    vector<Matrix2d> S(1);
    vector<Matrix2d> sigma(1);

    ConstitutiveModel.CauchyStress(false, mu0, lambda0, hardeningCoeff, sigma, JElastic, JPlastic,
                                    elasticDeformationGradient, R, S);

    // isApprox against an exact Zero() requires bit-exact equality (relative
    // tolerance degenerates when one side has zero norm), which is too
    // strict given the SVD inside ComputePolarDecomposition can leave R off
    // from I by floating-point noise. Compare the norm against a margin instead.
    REQUIRE(sigma[0].norm() == Catch::Approx(0.0).margin(1e-10));
}
