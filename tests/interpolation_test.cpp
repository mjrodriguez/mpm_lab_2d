#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/INTERPOLATION.h"

// Reference values below come directly from the closed-form cubic B-spline
// in Stomakhin et al. 2013, "A Material Point Method for Snow Simulation"
// (ACM TOG 32,4, Article 102), Section 4:
//
//   N(x) = { (1/2)|x|^3 - x^2 + 2/3,        0 <= |x| < 1
//          { -(1/6)|x|^3 + x^2 - 2|x| + 4/3, 1 <= |x| < 2
//          { 0,                              otherwise

TEST_CASE("INTERPOLATION::N matches the closed-form cubic B-spline", "[interpolation]") {
    INTERPOLATION Interpolation;

    REQUIRE(Interpolation.N(0.0) == Catch::Approx(2.0 / 3.0));
    REQUIRE(Interpolation.N(0.5) == Catch::Approx(0.5 * 0.125 - 0.25 + 2.0 / 3.0));
    REQUIRE(Interpolation.N(1.0) == Catch::Approx(1.0 / 6.0));
    REQUIRE(Interpolation.N(1.5) == Catch::Approx(-1.0 / 6.0 * 3.375 + 2.25 - 3.0 + 4.0 / 3.0));
    REQUIRE(Interpolation.N(2.0) == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(Interpolation.N(2.5) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("INTERPOLATION::N is even (symmetric about x = 0)", "[interpolation]") {
    INTERPOLATION Interpolation;

    for (double x : {0.2, 0.5, 0.9, 1.0, 1.3, 1.9, 2.5}) {
        REQUIRE(Interpolation.N(x) == Catch::Approx(Interpolation.N(-x)));
    }
}

TEST_CASE("INTERPOLATION::GradientN is the analytic derivative of N", "[interpolation]") {
    INTERPOLATION Interpolation;
    const double eps = 1e-6;

    for (double x : {0.3, 0.7, 1.3, 1.7, -0.4, -1.4}) {
        double numericalGradient = (Interpolation.N(x + eps) - Interpolation.N(x - eps)) / (2.0 * eps);
        REQUIRE(Interpolation.GradientN(x) == Catch::Approx(numericalGradient).margin(1e-6));
    }
}

TEST_CASE("INTERPOLATION::GradientN is odd (N is even, so N' is odd)", "[interpolation]") {
    INTERPOLATION Interpolation;

    for (double x : {0.2, 0.5, 0.9, 1.3, 1.9}) {
        REQUIRE(Interpolation.GradientN(x) == Catch::Approx(-Interpolation.GradientN(-x)));
    }
    REQUIRE(Interpolation.GradientN(0.0) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("INTERPOLATION::Weight (2D) forms a partition of unity over its support", "[interpolation]") {
    INTERPOLATION Interpolation;
    const double h = 1.0;
    Vector2d particlePosition(2.3, -1.7);

    double total = 0.0;
    for (int i = -3; i <= 6; i++) {
        for (int j = -6; j <= 3; j++) {
            total += Interpolation.Weight(h, particlePosition, i, j);
        }
    }

    REQUIRE(total == Catch::Approx(1.0));
}

TEST_CASE("INTERPOLATION::GradientWeight (2D) matches a central-difference gradient of Weight", "[interpolation]") {
    INTERPOLATION Interpolation;
    const double h = 1.0;
    const double eps = 1e-6;
    const double i = 2.0;
    const double j = -1.0;
    Vector2d particlePosition(2.6, -0.8);

    Vector2d analyticalGradient = Interpolation.GradientWeight(h, particlePosition, i, j);

    Vector2d xPlus = particlePosition + Vector2d(eps, 0);
    Vector2d xMinus = particlePosition - Vector2d(eps, 0);
    double dWdx = (Interpolation.Weight(h, xPlus, i, j) - Interpolation.Weight(h, xMinus, i, j)) / (2.0 * eps);

    Vector2d yPlus = particlePosition + Vector2d(0, eps);
    Vector2d yMinus = particlePosition - Vector2d(0, eps);
    double dWdy = (Interpolation.Weight(h, yPlus, i, j) - Interpolation.Weight(h, yMinus, i, j)) / (2.0 * eps);

    REQUIRE(analyticalGradient.x() == Catch::Approx(dWdx).margin(1e-6));
    REQUIRE(analyticalGradient.y() == Catch::Approx(dWdy).margin(1e-6));
}
