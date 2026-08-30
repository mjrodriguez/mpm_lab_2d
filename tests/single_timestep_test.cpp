#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "include/ELASTOPLASTIC.h"
#include "include/GRID.h"
#include "include/INTERPOLATION.h"
#include "include/PARTICLES.h"

// Integration smoke test: one ParticleToGrid -> ComputeGridForces ->
// UpdateVelocityGrid pass, for a small "free fall" particle cube (uniform
// initial velocity, no obstacle/collision -- GRID_NO_OBSTACLE), should
// conserve total momentum.
//
// This holds today because freshly-initialized particles have
// elasticDeformationGradient = Identity, which ELASTOPLASTIC::CauchyStress
// maps to exactly zero stress (see elastoplastic_test.cpp), so
// ComputeGridForces contributes zero force and UpdateVelocityGrid is a
// no-op (newVelocity = velocity + dt*0/mass = velocity). Note main.cpp's
// actual simulation loop never calls Grid.AddGravityForce (it's commented
// out there), so there is no external force to include in this test either
// -- logged as a separate finding in scratch/REVIEW.md.
//
// This test does NOT catch (and isn't meant to catch) the separate,
// higher-priority bug also logged in scratch/REVIEW.md: ParticleToGrid and
// ComputeGridForces both accumulate onto a single fixed grid node per
// particle instead of distributing across the B-spline neighborhood. Total
// momentum still balances either way (nothing is lost, just misplaced), so
// a momentum check alone can't distinguish correct distribution from that
// bug -- it needs a dedicated multi-node test instead.
TEST_CASE("One ParticleToGrid -> ComputeGridForces -> UpdateVelocityGrid pass conserves momentum (free fall, no collision)", "[integration]") {
    GRID_NO_OBSTACLE Grid;
    Grid.SetDefaultGrid();

    ELASTOPLASTIC ConstitutiveModel;
    INTERPOLATION Interpolation;

    RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
    Particle.m_numberOfParticles = 0;
    Particle.m_area = 0;
    Particle.SetInitialDensity(400.0);              // paper's default snow density (Table 2)
    Particle.SetCube(3, 0.05, Vector2d(0.5, 0.5));   // 9 particles, well inside the grid domain
    Particle.InitializeParticles();                  // mass/volume/density/F=I/velocity=(0,-3) via
                                                       // RECTANGLE_FREEFALL_HAT_OBSTACLE::InitializeParticleVelocities

    double totalParticleMass = 0.0;
    Vector2d initialMomentum = Vector2d::Zero();
    for (int p = 0; p < Particle.GetNumberOfParticles(); p++) {
        totalParticleMass += Particle.mass[p];
        initialMomentum += Particle.mass[p] * Particle.velocity[p];
    }
    REQUIRE(totalParticleMass == Catch::Approx(1.0));       // density * area = 400 * 0.05^2
    REQUIRE(initialMomentum.isApprox(Vector2d(0, -3)));     // uniform velocity (0, -3) * mass 1.0

    Grid.ParticleToGrid(Particle.mass, Particle.position, Particle.velocity, Interpolation);
    Grid.NodesWithMass();

    Vector2d gridMomentumAfterP2G = Vector2d::Zero();
    for (const Vector2i& node : Grid.massList) {
        int idx = Grid.Index2D(node.x(), node.y());
        gridMomentumAfterP2G += Grid.mass[idx] * Grid.velocity[idx];
    }
    REQUIRE(gridMomentumAfterP2G.isApprox(initialMomentum, 1e-9));

    Particle.ComputeVolumeDensity(Grid, Interpolation);  // matches main.cpp's first-iteration call

    const bool usePlasticity = false;
    const double mu0 = 100.0;
    const double lambda0 = 200.0;
    const double hardeningCoeff = 10.0;
    Grid.ComputeGridForces(usePlasticity, mu0, lambda0, hardeningCoeff, Particle.JElastic, Particle.JPlastic,
                            Particle.volume, Particle.position, Particle.cauchyStress,
                            Particle.elasticDeformationGradient, Particle.R, Particle.S, ConstitutiveModel,
                            Interpolation);

    for (const Vector2i& node : Grid.massList) {
        int idx = Grid.Index2D(node.x(), node.y());
        REQUIRE(Grid.force[idx].isApprox(Vector2d::Zero()));  // F_E = I everywhere -> zero stress -> zero force
    }

    const double dt = 1e-4;
    Grid.UpdateVelocityGrid(dt);

    Vector2d gridMomentumAfterUpdate = Vector2d::Zero();
    for (const Vector2i& node : Grid.massList) {
        int idx = Grid.Index2D(node.x(), node.y());
        gridMomentumAfterUpdate += Grid.mass[idx] * Grid.newVelocity[idx];
    }
    REQUIRE(gridMomentumAfterUpdate.isApprox(initialMomentum, 1e-9));
}

// Expected-failing (Catch2 "[!shouldfail]") regression test tracking the
// ParticleToGrid mass-distribution bug logged in scratch/REVIEW.md. Asserts
// the CORRECT cubic B-spline distribution: for a particle placed exactly
// halfway between grid nodes (a clean, symmetric case), mass should split
// evenly across the 4 nearest nodes, not land entirely on one. This
// currently fails because ParticleToGrid accumulates every node's weight
// into the fixed center index instead of Index2D(i,j) (GRID.cpp:265-266).
// Once that's fixed, this test should start passing -- at that point,
// remove the "[!shouldfail]" tag so it becomes a real regression guard
// instead of a tracked, known failure.
TEST_CASE("ParticleToGrid distributes mass across the B-spline neighborhood, not onto a single node", "[integration][!shouldfail]") {
    GRID_NO_OBSTACLE Grid;
    Grid.SetDefaultGrid();
    INTERPOLATION Interpolation;

    const double h = Grid.GetGridSpacing();
    const double particleMass = 2.0;
    const Vector2d position(0.5, 0.5);  // exactly halfway between two grid nodes in both x and y

    vector<double> massParticle = {particleMass};
    vector<Vector2d> positionParticle = {position};
    vector<Vector2d> velocityParticle = {Vector2d::Zero()};

    Grid.ParticleToGrid(massParticle, positionParticle, velocityParticle, Interpolation);

    int ig = static_cast<int>(std::floor((position.x() - Grid.Getxmin()) / h));
    int jg = static_cast<int>(std::floor((position.y() - Grid.Getymin()) / h));

    // By symmetry (particle exactly halfway between nodes), all four
    // nearest nodes get an equal, partial share of the particle's mass.
    double expectedWeightPerNode = Interpolation.Weight(h, position, ig, jg);
    REQUIRE(expectedWeightPerNode == Catch::Approx(0.229601).margin(1e-5));

    double expectedMassPerNode = particleMass * expectedWeightPerNode;

    REQUIRE(Grid.mass[Grid.Index2D(ig, jg)] == Catch::Approx(expectedMassPerNode));
    REQUIRE(Grid.mass[Grid.Index2D(ig + 1, jg)] == Catch::Approx(expectedMassPerNode));
    REQUIRE(Grid.mass[Grid.Index2D(ig, jg + 1)] == Catch::Approx(expectedMassPerNode));
    REQUIRE(Grid.mass[Grid.Index2D(ig + 1, jg + 1)] == Catch::Approx(expectedMassPerNode));
}
