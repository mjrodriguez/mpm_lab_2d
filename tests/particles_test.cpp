#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/PARTICLES.h"

// PARTICLES is abstract (pure virtual SetDefaultParticles/
// InitializeParticleVelocities/ParticleCollision), so tests instantiate a
// concrete subclass. SetCube/SetRectangle are defined on the base class and
// don't depend on which subclass is used -- RECTANGLE_FREEFALL_HAT_OBSTACLE
// is used here because it's the one main.cpp actually runs by default.
// (CUBE_FREEFALL, despite existing as a class declaration, has no method
// definitions anywhere in PARTICLES.cpp -- it's dead/unimplemented and
// fails to link if instantiated; found this while writing this test.)
//
// m_numberOfParticles and m_area are accumulated with += (so a simulation
// can compose multiple shapes, e.g. CUBE_TO_CUBE_FREEFALL), and are never
// initialized by a constructor -- every SetDefaultParticles() override
// zeroes them manually before calling SetCube/SetRectangle. Tests do the
// same to avoid reading uninitialized memory.

TEST_CASE("PARTICLES::SetCube places particlesPerDirection^2 particles on a square grid", "[particles]") {
    RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
    Particle.m_numberOfParticles = 0;
    Particle.m_area = 0;

    const double particlesPerDirection = 3;
    const double cubeLength = 1.2;
    const Vector2d anchor(0.1, 0.2);
    const double spacing = cubeLength / particlesPerDirection;  // 0.4

    Particle.SetCube(particlesPerDirection, cubeLength, anchor);

    REQUIRE(Particle.position.size() == 9);
    REQUIRE(Particle.GetNumberOfParticles() == Catch::Approx(9.0));
    REQUIRE(Particle.m_particleGridSpacing == Catch::Approx(spacing));
    REQUIRE(Particle.m_area == Catch::Approx(cubeLength * cubeLength));

    // position is filled in (px, py) order, index = px * N + py
    REQUIRE(Particle.position[0].isApprox(anchor));                                   // px=0, py=0
    REQUIRE(Particle.position[1].isApprox(anchor + Vector2d(0, spacing)));             // px=0, py=1
    REQUIRE(Particle.position[3].isApprox(anchor + Vector2d(spacing, 0)));             // px=1, py=0
    REQUIRE(Particle.position[8].isApprox(anchor + Vector2d(2 * spacing, 2 * spacing))); // px=2, py=2 (far corner)
}

TEST_CASE("PARTICLES::SetRectangle places particlesInX * particlesInY particles on a rectangular grid", "[particles]") {
    RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
    Particle.m_numberOfParticles = 0;
    Particle.m_area = 0;

    const double particlesInX = 2;
    const double particlesInY = 4;
    const double xLength = 1.0;
    const double yLength = 2.0;
    const Vector2d anchor(0.0, 0.0);
    const double xSpacing = xLength / particlesInX;  // 0.5
    const double ySpacing = yLength / particlesInY;  // 0.5

    Particle.SetRectangle(particlesInX, particlesInY, xLength, yLength, anchor);

    REQUIRE(Particle.position.size() == 8);
    REQUIRE(Particle.GetNumberOfParticles() == Catch::Approx(8.0));
    REQUIRE(Particle.m_area == Catch::Approx(xLength * yLength));

    // position is filled in (px, py) order, index = px * particlesInY + py
    REQUIRE(Particle.position[0].isApprox(anchor));                                       // px=0, py=0
    REQUIRE(Particle.position[1].isApprox(anchor + Vector2d(0, ySpacing)));                // px=0, py=1
    REQUIRE(Particle.position[4].isApprox(anchor + Vector2d(xSpacing, 0)));                // px=1, py=0
    REQUIRE(Particle.position[7].isApprox(anchor + Vector2d(xSpacing, 3 * ySpacing)));     // px=1, py=3 (far corner)
}

TEST_CASE("PARTICLES::SetCube/SetRectangle accumulate onto existing particle count and area", "[particles]") {
    RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
    Particle.m_numberOfParticles = 0;
    Particle.m_area = 0;

    Particle.SetCube(3, 1.2, Vector2d(0.1, 0.2));       // 9 particles, area 1.44
    Particle.SetRectangle(2, 4, 1.0, 2.0, Vector2d(0.0, 0.0));  // 8 particles, area 2.0

    REQUIRE(Particle.position.size() == 17);
    REQUIRE(Particle.GetNumberOfParticles() == Catch::Approx(17.0));
    REQUIRE(Particle.m_area == Catch::Approx(1.2 * 1.2 + 2.0));
}

// mass[p] = (density * area) / numberOfParticles for every particle, so
// total mass = density * area regardless of how many particles it's split
// across -- mass conservation. ("area" is this 2D simulation's stand-in for
// what a 3D version would call volume, per the TESTING_TODO wording.)
TEST_CASE("PARTICLES::InitializeParticleMass splits density * area evenly and conserves total mass", "[particles]") {
    const double density = 400.0;  // paper's default snow density (Table 2)

    SECTION("cube") {
        RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
        Particle.m_numberOfParticles = 0;
        Particle.m_area = 0;
        Particle.SetInitialDensity(density);
        Particle.SetCube(3, 1.2, Vector2d(0.1, 0.2));  // 9 particles, area 1.44

        Particle.InitializeParticleMass();

        const double totalMass = density * Particle.m_area;  // 576.0
        REQUIRE(Particle.mass.size() == 9);
        for (double m : Particle.mass) {
            REQUIRE(m == Catch::Approx(totalMass / 9.0));
        }

        double summedMass = 0.0;
        for (double m : Particle.mass) {
            summedMass += m;
        }
        REQUIRE(summedMass == Catch::Approx(totalMass));
    }

    SECTION("rectangle") {
        RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
        Particle.m_numberOfParticles = 0;
        Particle.m_area = 0;
        Particle.SetInitialDensity(density);
        Particle.SetRectangle(2, 4, 1.0, 2.0, Vector2d(0.0, 0.0));  // 8 particles, area 2.0

        Particle.InitializeParticleMass();

        const double totalMass = density * Particle.m_area;  // 800.0
        REQUIRE(Particle.mass.size() == 8);
        for (double m : Particle.mass) {
            REQUIRE(m == Catch::Approx(totalMass / 8.0));
        }

        double summedMass = 0.0;
        for (double m : Particle.mass) {
            summedMass += m;
        }
        REQUIRE(summedMass == Catch::Approx(totalMass));
    }
}
