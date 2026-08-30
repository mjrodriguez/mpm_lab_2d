#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/GRID.h"

// GRID has a pure virtual GridCollision, so tests instantiate the simplest
// concrete subclass. Index2D/Clamp don't depend on collision geometry.

TEST_CASE("GRID::Index2D maps (i, j) to a flat row-major index using m_ny", "[grid]") {
    GRID_NO_OBSTACLE Grid;
    Grid.m_ny = 5;

    REQUIRE(Grid.Index2D(0, 0) == 0);
    REQUIRE(Grid.Index2D(0, 4) == 4);
    REQUIRE(Grid.Index2D(1, 2) == 7);
    REQUIRE(Grid.Index2D(2, 3) == 13);
}

TEST_CASE("GRID::Clamp restricts principal stretches to [1 - criticalCompression, 1 + criticalStretch]", "[grid]") {
    GRID_NO_OBSTACLE Grid;
    const double criticalCompression = 0.2;  // lower bound: 1 - 0.2 = 0.8
    const double criticalStretch = 0.1;      // upper bound: 1 + 0.1 = 1.1

    SECTION("values outside the range are clamped, off-diagonal entries pass through") {
        Matrix2d principalValues;
        principalValues << 0.5, 0.3,
                            0.3, 1.5;

        Matrix2d clamped = Grid.Clamp(principalValues, criticalStretch, criticalCompression);

        REQUIRE(clamped(0, 0) == Catch::Approx(0.8));
        REQUIRE(clamped(1, 1) == Catch::Approx(1.1));
        REQUIRE(clamped(0, 1) == Catch::Approx(0.3));
        REQUIRE(clamped(1, 0) == Catch::Approx(0.3));

        // Clamp takes principalValues by non-const reference and mutates it
        // in place, in addition to returning a (now-identical) copy.
        REQUIRE(principalValues(0, 0) == Catch::Approx(0.8));
        REQUIRE(principalValues(1, 1) == Catch::Approx(1.1));
    }

    SECTION("values already inside the range are left unchanged") {
        Matrix2d principalValues;
        principalValues << 0.95, 0.0,
                            0.0, 1.05;

        Matrix2d clamped = Grid.Clamp(principalValues, criticalStretch, criticalCompression);

        REQUIRE(clamped(0, 0) == Catch::Approx(0.95));
        REQUIRE(clamped(1, 1) == Catch::Approx(1.05));
    }

    SECTION("values exactly on the boundary are left unchanged") {
        Matrix2d principalValues;
        principalValues << 0.8, 0.0,
                            0.0, 1.1;

        Matrix2d clamped = Grid.Clamp(principalValues, criticalStretch, criticalCompression);

        REQUIRE(clamped(0, 0) == Catch::Approx(0.8));
        REQUIRE(clamped(1, 1) == Catch::Approx(1.1));
    }
}
