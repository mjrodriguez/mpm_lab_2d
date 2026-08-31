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

// MakeParticles: CUBE_FREEFALL is deliberately NOT tested here (and not a
// valid name) -- it's declared in PARTICLES.h but has no method definitions
// in PARTICLES.cpp, so referencing `new CUBE_FREEFALL()` anywhere fails to
// link. See scratch/REVIEW.md.
TEST_CASE("MakeParticles returns the correct concrete type for each valid name", "[factory]") {
    SECTION("CUBE_TO_CUBE_COLLISION") {
        unique_ptr<PARTICLES> particles = MakeParticles("CUBE_TO_CUBE_COLLISION");
        REQUIRE(dynamic_cast<CUBE_TO_CUBE_COLLISION*>(particles.get()) != nullptr);
    }

    SECTION("RECTANGLE_FREEFALL_CYLINDER") {
        unique_ptr<PARTICLES> particles = MakeParticles("RECTANGLE_FREEFALL_CYLINDER");
        REQUIRE(dynamic_cast<RECTANGLE_FREEFALL_CYLINDER*>(particles.get()) != nullptr);
    }

    SECTION("RECTANGLE_FREEFALL_HAT_OBSTACLE") {
        unique_ptr<PARTICLES> particles = MakeParticles("RECTANGLE_FREEFALL_HAT_OBSTACLE");
        REQUIRE(dynamic_cast<RECTANGLE_FREEFALL_HAT_OBSTACLE*>(particles.get()) != nullptr);
    }

    SECTION("CUBE_CRASH_PADDED_GROUND") {
        unique_ptr<PARTICLES> particles = MakeParticles("CUBE_CRASH_PADDED_GROUND");
        REQUIRE(dynamic_cast<CUBE_CRASH_PADDED_GROUND*>(particles.get()) != nullptr);
    }

    SECTION("CUBE_CRASH_WALL") {
        unique_ptr<PARTICLES> particles = MakeParticles("CUBE_CRASH_WALL");
        REQUIRE(dynamic_cast<CUBE_CRASH_WALL*>(particles.get()) != nullptr);
    }

    SECTION("CUBE_TO_CUBE_FREEFALL") {
        unique_ptr<PARTICLES> particles = MakeParticles("CUBE_TO_CUBE_FREEFALL");
        REQUIRE(dynamic_cast<CUBE_TO_CUBE_FREEFALL*>(particles.get()) != nullptr);
    }
}

TEST_CASE("MakeParticles throws a clear error for an unknown name", "[factory]") {
    REQUIRE_THROWS_AS(MakeParticles("NOT_A_REAL_PARTICLES"), runtime_error);
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("NOT_A_REAL_PARTICLES"));

    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("CUBE_TO_CUBE_COLLISION"));
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("RECTANGLE_FREEFALL_CYLINDER"));
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("RECTANGLE_FREEFALL_HAT_OBSTACLE"));
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("CUBE_CRASH_PADDED_GROUND"));
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("CUBE_CRASH_WALL"));
    REQUIRE_THROWS_WITH(MakeParticles("NOT_A_REAL_PARTICLES"), Catch::Matchers::ContainsSubstring("CUBE_TO_CUBE_FREEFALL"));
}

TEST_CASE("MakeParticles throws on CUBE_FREEFALL (deliberately excluded, dead class)", "[factory]") {
    REQUIRE_THROWS_AS(MakeParticles("CUBE_FREEFALL"), runtime_error);
}

TEST_CASE("MakeParticles throws on an empty name", "[factory]") {
    REQUIRE_THROWS_AS(MakeParticles(""), runtime_error);
}

TEST_CASE("MakeSimulationParameters returns the correct concrete type for each valid name", "[factory]") {
    SECTION("DEFAULT_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("DEFAULT_PARAMETERS");
        REQUIRE(dynamic_cast<DEFAULT_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("DEFAULT_IMPLICIT_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("DEFAULT_IMPLICIT_PARAMETERS");
        REQUIRE(dynamic_cast<DEFAULT_IMPLICIT_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("LOWER_YOUNGS_MODULUS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_YOUNGS_MODULUS");
        REQUIRE(dynamic_cast<LOWER_YOUNGS_MODULUS*>(params.get()) != nullptr);
    }

    SECTION("LOWER_YOUNGS_MODULUS_IMPLICIT") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_YOUNGS_MODULUS_IMPLICIT");
        REQUIRE(dynamic_cast<LOWER_YOUNGS_MODULUS_IMPLICIT*>(params.get()) != nullptr);
    }

    SECTION("LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS");
        REQUIRE(dynamic_cast<LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("LOWER_HARDENING") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_HARDENING");
        REQUIRE(dynamic_cast<LOWER_HARDENING*>(params.get()) != nullptr);
    }

    SECTION("LOWER_CRITICAL_COMPRESSION_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_CRITICAL_COMPRESSION_PARAMETERS");
        REQUIRE(dynamic_cast<LOWER_CRITICAL_COMPRESSION_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("LOWER_CRITICAL_STRETCH_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("LOWER_CRITICAL_STRETCH_PARAMETERS");
        REQUIRE(dynamic_cast<LOWER_CRITICAL_STRETCH_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("HYPERELASTICITY") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("HYPERELASTICITY");
        REQUIRE(dynamic_cast<HYPERELASTICITY*>(params.get()) != nullptr);
    }

    SECTION("SCALED_DEFAULT_PARAMETERS") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("SCALED_DEFAULT_PARAMETERS");
        REQUIRE(dynamic_cast<SCALED_DEFAULT_PARAMETERS*>(params.get()) != nullptr);
    }

    SECTION("SCALED_HYPERELASTICITY") {
        unique_ptr<SIMULATION_PARAMETERS> params = MakeSimulationParameters("SCALED_HYPERELASTICITY");
        REQUIRE(dynamic_cast<SCALED_HYPERELASTICITY*>(params.get()) != nullptr);
    }
}

TEST_CASE("MakeSimulationParameters throws a clear error for an unknown name", "[factory]") {
    REQUIRE_THROWS_AS(MakeSimulationParameters("NOT_REAL_PARAMETERS"), runtime_error);
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("NOT_REAL_PARAMETERS"));

    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("DEFAULT_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("DEFAULT_IMPLICIT_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_YOUNGS_MODULUS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_YOUNGS_MODULUS_IMPLICIT"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_HARDENING"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_CRITICAL_COMPRESSION_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("LOWER_CRITICAL_STRETCH_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("HYPERELASTICITY"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("SCALED_DEFAULT_PARAMETERS"));
    REQUIRE_THROWS_WITH(MakeSimulationParameters("NOT_REAL_PARAMETERS"), Catch::Matchers::ContainsSubstring("SCALED_HYPERELASTICITY"));
}

TEST_CASE("MakeSimulationParameters throws on an empty name", "[factory]") {
    REQUIRE_THROWS_AS(MakeSimulationParameters(""), runtime_error);
}
