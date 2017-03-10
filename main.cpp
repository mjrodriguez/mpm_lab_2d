//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <iostream>
#include <vector>


#include "include/ELASTOPLASTIC.h"
#include "include/FILE_IO.h"
#include "include/GRID.h"
#include "include/INTERPOLATION.h"
#include "include/PARTICLES.h"
#include "include/SIMULATION_PARAMETERS.h"
#include "include/TOOLS.h"

using namespace std;
using namespace Eigen;

int main() {

    ELASTOPLASTIC ConstitutiveModel;
    FILE_IO FileIO;
    GRID Grid;
    INTERPOLATION Interpolation;
    PARTICLES Particle;
    DEFAULT_PARAMETERS SimulationParams;
    // HYPERELASTICITY SimulationParams;
	// SCALED_HYPERELASTICITY SimulationParams;
	// SCALED_DEFAULT_PARAMETERS SimulationParams;
    // LOWER_CRITICAL_COMPRESSION_PARAMETERS SimulationParams;
    TOOLS Tools;


    string directory = "/home/rodriguez/Documents/mpm_lab_node_4/output/";
    // string directory = "/home/rodriguez/Documents/mpm_lab_node_9/output/";
    // string directory = "/Users/Martin/Documents/College/Research/Professor Blomgren/Thesis/mpm_lab";

    SimulationParams.SetDefaultParameters();
    SimulationParams.SetDt(5e-5);
    Particle.SetInitialVelocity(Vector2d(0.0,-1.0));
    Particle.SetDefaultParticles();
    Grid.SetDefaultGrid();

    cout << "Run this simulation? " << SimulationParams.GetSimulationName() << endl;
    cout << "Save to this directory? " << directory << endl;
    cin.get();

    ///////////////////////////////////////
    // SAVE SIMULATION PARAMETERS
    ////////////////////////////////////////
    FileIO.WriteSimulationParameters(directory, Particle, Grid, SimulationParams);


    double currentTime = 0;
    int iterationCounter = 0;
    int frameNumber = 0;
    vector<double> timePerIteration;
    vector<double> timeStep;
    timeStep.push_back(SimulationParams.GetDt());
	
    // WRITING INITIAL CONDITIONS
    FileIO.WriteOutput(directory, SimulationParams.GetSimulationName(), frameNumber, Particle.position, Particle.velocity, Particle.cauchyStress, Particle.deformationGradient, Particle.elasticDeformationGradient, Particle.plasticDeformationGradient);
	
    clock_t begin_sim = clock();

    while (currentTime <= SimulationParams.GetFinalTime()  ){
        clock_t begin = clock();

        cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << "Simulation = " << SimulationParams.GetSimulationName() << " | Num of Particles = " << Particle.GetNumberOfParticles() << " | Current Time Step Size: " << SimulationParams.GetDt() << " | " << "Current Time: " << currentTime << " | Current Iteration: " << iterationCounter << " | Current Frame: " << frameNumber <<  endl;
        cout << "Saving Directory = " << directory << endl;

        iterationCounter += 1;

        Grid.ComputeGridNodeWeights(Particle.position, Interpolation);
        Grid.ParticleToGrid( Particle.mass, Particle.position, Particle.velocity, Interpolation);
        Grid.NodesWithMass();

        cout << "Mass on Particles = " << Tools.Sum(Particle.mass) << endl;
        cout << "Mass on Grid = " << Tools.Sum(Grid.mass) << endl;

        if (currentTime == 0){
            Particle.ComputeVolumeDensity(Grid.GetN(), Grid.GetGridSpacing(), Grid.mass, Grid.massList, Tools, Interpolation);
        }


        Grid.ComputeGridForces(SimulationParams.usePlasticity, SimulationParams.GetMu(), SimulationParams.GetLambda(), SimulationParams.GetHardeningCoeff(), Particle.JElastic, Particle.JPlastic, Particle.volume, Particle.position, Particle.cauchyStress, Particle.elasticDeformationGradient, Particle.R, Particle.S, ConstitutiveModel, Interpolation);

        cout << "Max Force from Internal Stress = " << Tools.MaxNormValue(Grid.force) << endl;
        cout << "Min Force from Internal Stress = " << Tools.MinNormValue(Grid.force) << endl;

        Grid.AddGravityForce(SimulationParams.gravity);
        Grid.UpdateVelocityGrid(SimulationParams.GetDt());


        Grid.GridCollision(SimulationParams.GetFrictionCoeff());
        Grid.UpdateDeformationGradient(SimulationParams.usePlasticity, SimulationParams.GetDt(), SimulationParams.GetCriticalStretch(),
                                       SimulationParams.GetCriticalCompression(), Particle.JElastic, Particle.JPlastic,
                                       Particle.position, Particle.elasticDeformationGradient, Particle.plasticDeformationGradient, Particle.deformationGradient, Interpolation);


        Particle.ComputePicFlipVelocities(Grid.GetN(), Grid.GetGridSpacing(), Grid.massList, Grid.nodeWeights, Grid.velocity, Grid.newVelocity, Tools, Interpolation);
        cout << "Max PIC Velocity = " << Tools.MaxNormValue(Particle.velocityPIC) << endl;
        cout << "Max FLIP Velocity = " << Tools. MaxNormValue(Particle.velocityFLIP) << endl;
		
		Particle.UpdateParticleVelocity(SimulationParams.GetDt());
        Particle.ParticleCollision(SimulationParams.GetFrictionCoeff());
		Particle.UpdateParticlePosition(SimulationParams.GetDt());

        cout << "Max Velocity = " << Tools.MaxNormValue(Particle.velocity) << endl;
        cout << "Cube 1 = " << Particle.position[Particle.GetNumberOfParticles()/4].transpose() << endl;
		cout << "Cube 2 = " << Particle.position[3*Particle.GetNumberOfParticles()/4].transpose() << endl;
		
        if ( iterationCounter % 500 == 0){
        //if (1 == 0){
            frameNumber += 1;

            // WRITING OUTPUT
            FileIO.WriteOutput(directory, SimulationParams.GetSimulationName(), frameNumber, Particle.position, Particle.velocity, Particle.cauchyStress, Particle.deformationGradient, Particle.elasticDeformationGradient, Particle.plasticDeformationGradient);



//			// COMPUTING DISTANCE TRAVELED OF PARTICLE
//			distanceTraveled = Particle.position.back() - oldPosition;
//			oldPosition = Particle.position.back();
//			cout << "Distance Traveled.norm() = " << distanceTraveled.norm() << endl;

            // cin.get();
        }

        Grid.ClearEulerianFields();

        currentTime += SimulationParams.GetDt();

        //SimulationParams.SetDt(Grid.GetGridSpacing(), Particle.velocity );
        timeStep.push_back(SimulationParams.GetDt());

        clock_t end = clock();
        timePerIteration.push_back( double(end - begin)/CLOCKS_PER_SEC );
        double iterationTimeDisplay = double(end - begin)/CLOCKS_PER_SEC;
        cout << "Time per Iteration = " << iterationTimeDisplay << " seconds" << endl;


    }

    clock_t end_sim = clock();
    double totalTime = double(end_sim - begin_sim)/CLOCKS_PER_SEC;
    cout << "Time for Simulation = " << totalTime << " seconds." << endl;

    ///////////////////////////////////////////////////////////////////////////////
    // WRITING TIMESTEP AND TIME PER ITERATION DATA
    ///////////////////////////////////////////////////////////////////////////////

    FileIO.WriteTimeStep(directory, SimulationParams.GetSimulationName(), timeStep);
    FileIO.WriteTimePerIteration(directory, SimulationParams.GetSimulationName(), timePerIteration);


    // std::cout << "Hello, World!" << std::endl;
    return 0;








}