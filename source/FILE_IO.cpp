//
// Created by Martin Rodriguez Cruz on 3/7/17.
//


//
// Created by Martin Rodriguez Cruz on 9/18/16.
//


#include <vector>
#include "../Eigen/Dense"
#include <fstream>
#include <iostream>

#include "../include/FILE_IO.h"
#include "../include/SIMULATION_PARAMETERS.h"
#include "../include/PARTICLES.h"
#include "../include/GRID.h"

using namespace std;
using namespace Eigen;


////////////////////////////////////////////////////////////////////////
//*********************************************************************
// WRITING SIMULATION PARAMETERS TO TEXT FILE
// Author: Martin Rodriguez Cruz
// Date: February 21, 2017
// Input: 	directory --> saving directory, should be string
// 			Particle --> Pass particle class
//			SimulationParameters --> Pass parameter class
//*********************************************************************
///////////////////////////////////////////////////////////////////////

void FILE_IO::WriteSimulationParameters(string directory, PARTICLES& Particles, GRID& Grid, SIMULATION_PARAMETERS& SimulationParameters){
	
    string name = SimulationParameters.GetSimulationName() + string("_") + string("simulation_parameters");

    string filename(directory + name + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());

    fileout << "Name of Simulation = " << SimulationParameters.GetSimulationName() << endl;
    fileout << "Initial Velocity (Particles) = " << Particles.GetInitialVelocity().transpose() << endl;
	fileout << "Number of Particles = " << Particles.GetNumberOfParticles() << endl;
    fileout << "Grid size = " << Grid.GetN() << endl;
    fileout << "Grid Spacing = " << Grid.GetGridSpacing() << endl;
	fileout << "Timestep Length = " << SimulationParameters.GetDt() << endl;
    fileout << "Critical Compression = " << SimulationParameters.GetCriticalCompression() << endl;
    fileout << "Critical Stretch = " << SimulationParameters.GetCriticalStretch() << endl;
    fileout << "Hardening Coefficient = " << SimulationParameters.GetHardeningCoeff() << endl;
    fileout << "Initial Density = " << SimulationParameters.GetInitialDensity() << endl;
    fileout << "Young's Modulus = " << SimulationParameters.GetYoungsMod() << endl;
    fileout << "Poisson Ratio = " << SimulationParameters.GetPoisson() << endl;
    fileout << "Lambda0 = " << SimulationParameters.GetLambda() << endl;
    fileout << "Mu0 = " << SimulationParameters.GetMu() << endl;
    fileout << "Using Plasticity = " << SimulationParameters.usePlasticity << endl;
    fileout << "Interpolation alpha = " << SimulationParameters.GetAlpha() << endl;

    fileout.close();

}







/*****************************************************************************************
* WRITING OUTPUT TO TEXT FILES
****************************************************************************************/

void FILE_IO::WriteParticlePositions(string directory, string simulationName, vector<Vector2d> &particlePosition,
                                     int frameNumber) {

    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("particle_position");

    string filename(directory + string("/particle_position/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < particlePosition.size(); p++) {

        for (int i = 0; i < 2; i++) {
            fileout << particlePosition[p](i) << endl;
        }

    }

    fileout.close();
}

void FILE_IO::WriteParticleVelocities(string directory, string simulationName, vector<Vector2d> &particleVelocity,
                                      int frameNumber) {

    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("particle_velocity");

    string filename(directory + string("/particle_velocity/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < particleVelocity.size(); p++) {

        for (int i = 0; i < 2; i++) {
            fileout << particleVelocity[p](i) << endl;
        }

    }

    fileout.close();

}


void FILE_IO::WriteCauchyStress(string directory, string simulationName,  vector<Matrix2d> &cauchyStress, int frameNumber) {
    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("cauchy_stress");

    string filename(directory + string("/cauchy_stress/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < cauchyStress.size(); p++) {

        fileout << cauchyStress[p] << endl;

    }
    fileout.close();
}


void FILE_IO::WritedPsidF(string directory, string simulationName, vector<Matrix2d> &dpsi_df, int frameNumber) {
    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("dpsi_df");

    string filename(directory + string("/dPsi_dF/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < dpsi_df.size(); p++) {

        fileout << dpsi_df[p] << endl;

    }
}


void FILE_IO::WriteDeformationGradient(string directory, string simulationName, vector<Matrix2d> &deformationGradient,
                                       int frameNumber) {

    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("deformation_gradient");

    string filename(directory + string("/deformation_gradient/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < deformationGradient.size(); p++) {

        fileout << deformationGradient[p] << endl;

    }
    fileout.close();
}


void FILE_IO::WriteElasticDeformationGradient(string directory, string simulationName,
                                              vector<Matrix2d> &elasticDeformationGradient, int frameNumber) {
    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("elastic_deformation_gradient");

    string filename(directory + string("/elastic_deformation_gradient/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < elasticDeformationGradient.size(); p++) {

        fileout << elasticDeformationGradient[p] << endl;

    }
    fileout.close();
}


void FILE_IO::WritePlasticDeformationGradient(string directory, string simulationName,
                                              vector<Matrix2d> &plasticDeformationGradient, int frameNumber) {
    ostringstream iter;
    iter << "_" << frameNumber;

    string name = simulationName + string("_") + string("plastic_deformation_gradient");

    string filename(directory + string("/plastic_deformation_gradient/") + name + iter.str() + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < plasticDeformationGradient.size(); p++) {

        fileout << plasticDeformationGradient[p] << endl;

    }

    fileout.close();
}



void FILE_IO::WriteTimePerIteration(string directory, string simulationName, vector<double>& timePerIteration) {

    string name = simulationName + string("_") + string("time_per_iteration");

    string filename(directory + string("/") + name + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < timePerIteration.size(); p++) {
      fileout << timePerIteration[p] << endl;
    }

    fileout.close();

}

void FILE_IO::WriteTimeStep(string directory, string simulationName, vector<double>& timeStep) {

    string name = simulationName + string("_") + string("time_step");

    string filename(directory + string("/") + name + string(".txt"));
    ofstream fileout(filename.c_str());

    assert(fileout.is_open());


    for (int p = 0; p < timeStep.size(); p++) {
      fileout << timeStep[p] << endl;
    }

    fileout.close();
}






void FILE_IO::WriteOutput(string directory, string simulationName, int frameNumber, vector<Vector2d>& particlePosition, vector<Vector2d>& particleVelocity, vector<Matrix2d>& cauchyStress, vector<Matrix2d>& deformationGradient, vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& plasticDeformationGradient){
	
	WriteParticlePositions(directory, simulationName, particlePosition, frameNumber);
	WriteParticleVelocities(directory, simulationName, particleVelocity, frameNumber);
	WriteCauchyStress(directory, simulationName, cauchyStress, frameNumber);
	WriteDeformationGradient(directory, simulationName, deformationGradient, frameNumber);
	WriteElasticDeformationGradient(directory, simulationName, elasticDeformationGradient, frameNumber);
	WritePlasticDeformationGradient(directory, simulationName, plasticDeformationGradient, frameNumber);
	
}