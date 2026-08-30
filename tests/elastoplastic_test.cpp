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
