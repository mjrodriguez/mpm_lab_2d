#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "include/CONFIG.h"

// MakeGrid should return the specific concrete subclass named by the
// config, not just "any GRID" -- dynamic_cast to the expected type confirms
// that, rather than only checking the pointer is non-null.
TEST_CASE("MakeGrid returns the correct concrete type for each valid name", "[factory]") {
    SECTION("GRID_NO_OBSTACLE") {
        unique_ptr<GRID> grid = MakeGrid("GRID_NO_OBSTACLE");
        REQUIRE(dynamic_cast<GRID_NO_OBSTACLE*>(grid.get()) != nullptr);
    }

    SECTION("GRID_CYLINDER_OBSTACLE") {
        unique_ptr<GRID> grid = MakeGrid("GRID_CYLINDER_OBSTACLE");
        REQUIRE(dynamic_cast<GRID_CYLINDER_OBSTACLE*>(grid.get()) != nullptr);
    }

    SECTION("GRID_HAT_OBSTACLE") {
        unique_ptr<GRID> grid = MakeGrid("GRID_HAT_OBSTACLE");
        REQUIRE(dynamic_cast<GRID_HAT_OBSTACLE*>(grid.get()) != nullptr);
    }
}

TEST_CASE("MakeGrid throws a clear error for an unknown name", "[factory]") {
    REQUIRE_THROWS_AS(MakeGrid("NOT_A_REAL_GRID"), runtime_error);

    REQUIRE_THROWS_WITH(MakeGrid("NOT_A_REAL_GRID"), Catch::Matchers::ContainsSubstring("NOT_A_REAL_GRID"));

    // The error should list all three valid options, so a user hitting
    // this doesn't have to go dig through the source to find them.
    REQUIRE_THROWS_WITH(MakeGrid("NOT_A_REAL_GRID"), Catch::Matchers::ContainsSubstring("GRID_NO_OBSTACLE"));
    REQUIRE_THROWS_WITH(MakeGrid("NOT_A_REAL_GRID"), Catch::Matchers::ContainsSubstring("GRID_CYLINDER_OBSTACLE"));
    REQUIRE_THROWS_WITH(MakeGrid("NOT_A_REAL_GRID"), Catch::Matchers::ContainsSubstring("GRID_HAT_OBSTACLE"));
}

TEST_CASE("MakeGrid throws on an empty name", "[factory]") {
    REQUIRE_THROWS_AS(MakeGrid(""), runtime_error);
}
