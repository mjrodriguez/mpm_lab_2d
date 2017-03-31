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
	INTERPOLATION Interpolation;
    
	// GRID_NO_OBSTACLE Grid;
	GRID_HAT_OBSTACLE Grid;
	// GRID_CYLINDER_OBSTACLE Grid;
	
    // CUBE_TO_CUBE_FREEFALL Particle;
   	// CUBE_CRASH_WALL Particle;
	// CUBE_TO_CUBE_COLLISION Particle;
    // CUBE_CRASH_PADDED_GROUND Particle;
	RECTANGLE_FREEFALL_HAT_OBSTACLE Particle;
	// RECTANGLE_FREEFALL_CYLINDER Particle;
	
	DEFAULT_PARAMETERS SimulationParams;
    // HYPERELASTICITY SimulationParams;
    // LOWER_CRITICAL_COMPRESSION_PARAMETERS SimulationParams;
    // LOWER_CRITICAL_STRETCH_PARAMETERS SimulationParams;
    // LOWER_HARDENING SimulationParams;
    // LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS SimulationParams;
    // LOWER_YOUNGS_MODULUS SimulationParams;
    TOOLS Tools;


	string simulationNumber = "simulation_125";
    string directory = "/home/rodriguez/Documents/mpm_lab_node_5/output/";
    // string directory = "/home/rodriguez/Documents/mpm_lab_node_9/output/";
    // string directory = "/Users/Martin/Documents/College/Research/Professor Blomgren/Thesis/mpm_lab";

    SimulationParams.SetDefaultParameters();
	Particle.SetInitialDensity(SimulationParams.GetInitialDensity());
    Particle.SetDefaultParticles();
    Grid.SetDefaultGrid();
	double timeToFrame = SimulationParams.GetTimeToFrame();
	SimulationParams.SetDt(Grid.GetGridSpacing(), timeToFrame, Particle.velocity );
	
	
	SimulationParams.SetSimulationName(simulationNumber + "_" + SimulationParams.GetSimulationName());
		
    cout << "Run this simulation? " << SimulationParams.GetSimulationName() << "_"+Particle.GetParticleSimulationName() << endl;
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

        cout << "-------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << "Simulation = " << SimulationParams.GetSimulationName() << ", " << Particle.GetParticleSimulationName() << " | Num of Particles = " << Particle.GetNumberOfParticles() << " | Current Time Step Size: " << SimulationParams.GetDt() << endl; 
		cout << "Current Time: " << currentTime << " | Current Iteration: " << iterationCounter << " | Current Frame: " << frameNumber <<  endl;
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
		
		//////////////////////////////////////////////////////////////////
		// COMPUTE newVelocity from velocity
		/////////////////////////////////////////////////////////////////
        Grid.UpdateVelocityGrid(SimulationParams.GetDt());
		
		/////////////////////////////////////////////////////////////////////////
		// EVERYTHING IS UPDATED WITH newVelocity (not velocity)
		/////////////////////////////////////////////////////////////////////////

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
        cout << "First Particle  = " << Particle.position[0].transpose() << endl;
		// cout << "Cube 2 = " << Particle.position[3*Particle.GetNumberOfParticles()/4].transpose() << endl;
		
		timeToFrame -= SimulationParams.GetDt();
		
		if ( timeToFrame == 0){
			frameNumber += 1;

			// WRITING OUTPUT
			FileIO.WriteOutput(directory, SimulationParams.GetSimulationName(), frameNumber, Particle.position, Particle.velocity, Particle.cauchyStress, Particle.deformationGradient, Particle.elasticDeformationGradient, Particle.plasticDeformationGradient);
			timeToFrame = SimulationParams.GetTimeToFrame();
		}

        Grid.ClearEulerianFields();

        currentTime += SimulationParams.GetDt();
		
		
        SimulationParams.SetDt(Grid.GetGridSpacing(), timeToFrame, Particle.velocity );
        timeStep.push_back(SimulationParams.GetDt());

        clock_t end = clock();
        timePerIteration.push_back( double(end - begin)/CLOCKS_PER_SEC );
        double iterationTimeDisplay = double(end - begin)/CLOCKS_PER_SEC;
        cout << "Time per Iteration = " << iterationTimeDisplay << " seconds" << endl;
		cout << "Time to Frame = " << timeToFrame << " seconds" << endl;
		
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