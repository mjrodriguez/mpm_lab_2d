#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/GRID.h"
#include "include/PARTICLES.h"

// GradientLevelSet returns the *normalized* gradient of LevelSet (a unit
// surface normal), not the raw gradient -- so the central-difference
// estimate has to be normalized too before comparing. The GRID and
// PARTICLES classes below duplicate the same LevelSet/GradientLevelSet
// formulas (once per side of the collision handling), so both get tested.

namespace {

template <typename LevelSetFn>
Vector2d CentralDifferenceGradient(LevelSetFn levelSet, double x, double y, double eps = 1e-6) {
    double dfdx = (levelSet(x + eps, y) - levelSet(x - eps, y)) / (2.0 * eps);
    double dfdy = (levelSet(x, y + eps) - levelSet(x, y - eps)) / (2.0 * eps);
    Vector2d gradient(dfdx, dfdy);
    return gradient.normalized();
}

void RequireApproxEqual(const Vector2d& a, const Vector2d& b, double margin = 1e-5) {
    REQUIRE(a.x() == Catch::Approx(b.x()).margin(margin));
    REQUIRE(a.y() == Catch::Approx(b.y()).margin(margin));
}

}  // namespace

TEST_CASE("GRID_CYLINDER_OBSTACLE::GradientLevelSet matches the normalized central-difference gradient", "[levelset]") {
    GRID_CYLINDER_OBSTACLE Grid;
    const double x = 0.55;
    const double y = 0.35;  // away from the circle's center (0.5, 0.3), where the gradient is singular

    Vector2d numericalGradient = CentralDifferenceGradient(
        [&](double px, double py) { return Grid.LevelSet(px, py); }, x, y);

    RequireApproxEqual(Grid.GradientLevelSet(x, y), numericalGradient);
}

TEST_CASE("RECTANGLE_FREEFALL_CYLINDER::GradientLevelSet matches the normalized central-difference gradient", "[levelset]") {
    RECTANGLE_FREEFALL_CYLINDER Particle;
    const double x = 0.55;
    const double y = 0.35;

    Vector2d numericalGradient = CentralDifferenceGradient(
        [&](double px, double py) { return Particle.LevelSet(px, py); }, x, y);

    RequireApproxEqual(Particle.GradientLevelSet(x, y), numericalGradient);
}

TEST_CASE("GRID_HAT_OBSTACLE::GradientLevelSet matches the normalized central-difference gradient", "[levelset]") {
    GRID_HAT_OBSTACLE Grid;

    // Points are kept strictly inside the hat region (x in [0.4, 0.6], y in
    // [0.5, 0.6]) and away from the x = 0.5 ridge, where |x - 0.5| has a
    // kink and the central difference would straddle two different formulas.
    SECTION("left half of the hat (x < 0.5)") {
        const double x = 0.42;
        const double y = 0.55;

        Vector2d numericalGradient = CentralDifferenceGradient(
            [&](double px, double py) { return Grid.LevelSet(px, py); }, x, y);

        RequireApproxEqual(Grid.GradientLevelSet(x, y), numericalGradient);
    }

    SECTION("right half of the hat (x > 0.5)") {
        const double x = 0.58;
        const double y = 0.55;

        Vector2d numericalGradient = CentralDifferenceGradient(
            [&](double px, double py) { return Grid.LevelSet(px, py); }, x, y);

        RequireApproxEqual(Grid.GradientLevelSet(x, y), numericalGradient);
    }
}

TEST_CASE("RECTANGLE_FREEFALL_HAT_OBSTACLE::GradientLevelSet matches the normalized central-difference gradient", "[levelset]") {
    RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;

    SECTION("left half of the hat (x < 0.5)") {
        const double x = 0.42;
        const double y = 0.55;

        Vector2d numericalGradient = CentralDifferenceGradient(
            [&](double px, double py) { return Particle.LevelSet(px, py); }, x, y);

        RequireApproxEqual(Particle.GradientLevelSet(x, y), numericalGradient);
    }

    SECTION("right half of the hat (x > 0.5)") {
        const double x = 0.58;
        const double y = 0.55;

        Vector2d numericalGradient = CentralDifferenceGradient(
            [&](double px, double py) { return Particle.LevelSet(px, py); }, x, y);

        RequireApproxEqual(Particle.GradientLevelSet(x, y), numericalGradient);
    }
}
