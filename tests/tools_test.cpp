#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/TOOLS.h"

TEST_CASE("TOOLS::Sum totals a vector of doubles", "[tools]") {
    SECTION("empty vector sums to zero") {
        vector<double> input;
        REQUIRE(TOOLS::Sum(input) == 0.0);
    }

    SECTION("positive values") {
        vector<double> input = {1.0, 2.0, 3.0};
        REQUIRE(TOOLS::Sum(input) == Catch::Approx(6.0));
    }

    SECTION("mixed sign values") {
        vector<double> input = {1.5, -2.5, 3.0};
        REQUIRE(TOOLS::Sum(input) == Catch::Approx(2.0));
    }
}

TEST_CASE("TOOLS::Index2D maps (i, j) to a flat row-major index", "[tools]") {
    REQUIRE(TOOLS::Index2D(5, 0, 0) == 0);
    REQUIRE(TOOLS::Index2D(5, 1, 2) == 7);
    REQUIRE(TOOLS::Index2D(4, 3, 3) == 15);
}

TEST_CASE("TOOLS::Index3D maps (i, j, k) to a flat row-major index", "[tools]") {
    REQUIRE(TOOLS::Index3D(3, 0, 0, 0) == 0);
    REQUIRE(TOOLS::Index3D(3, 1, 1, 1) == 13);
    REQUIRE(TOOLS::Index3D(4, 2, 3, 1) == 45);
}

TEST_CASE("TOOLS::MaxNormValue/MinNormValue for Vector2d", "[tools]") {
    vector<Vector2d> V = {Vector2d(3, 4), Vector2d(0, 0), Vector2d(1, 1)};

    REQUIRE(TOOLS::MaxNormValue(V) == Catch::Approx(5.0));
    REQUIRE(TOOLS::MinNormValue(V) == Catch::Approx(0.0));
}

TEST_CASE("TOOLS::MaxNormValue/MinNormValue for Vector3d", "[tools]") {
    vector<Vector3d> V = {Vector3d(1, 2, 2), Vector3d(0, 0, 0), Vector3d(3, 4, 0)};

    REQUIRE(TOOLS::MaxNormValue(V) == Catch::Approx(5.0));
    REQUIRE(TOOLS::MinNormValue(V) == Catch::Approx(0.0));
}

TEST_CASE("TOOLS::MaxNormValue/MinNormValue for Matrix2d use the Frobenius norm", "[tools]") {
    Matrix2d zero = Matrix2d::Zero();
    Matrix2d identity = Matrix2d::Identity();
    Matrix2d threeFour;
    threeFour << 3, 4, 0, 0;

    vector<Matrix2d> V = {identity, zero, threeFour};

    REQUIRE(TOOLS::MaxNormValue(V) == Catch::Approx(5.0));
    REQUIRE(TOOLS::MinNormValue(V) == Catch::Approx(0.0));
}

TEST_CASE("TOOLS::MaxNormValue/MinNormValue for Matrix3d use the Frobenius norm", "[tools]") {
    Matrix3d zero = Matrix3d::Zero();
    Matrix3d identity = Matrix3d::Identity();
    Matrix3d diag345 = Vector3d(3, 4, 0).asDiagonal();

    vector<Matrix3d> V = {identity, zero, diag345};

    REQUIRE(TOOLS::MaxNormValue(V) == Catch::Approx(5.0));
    REQUIRE(TOOLS::MinNormValue(V) == Catch::Approx(0.0));
}
