//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#ifndef MPM_LAB_2D_FILE_IO_H
#define MPM_LAB_2D_FILE_IO_H

#include <vector>
#include "../Eigen/Dense"
#include "PARTICLES.h"
#include "SIMULATION_PARAMETERS.h"
#include "GRID.h"


using namespace std;
using namespace Eigen;


class FILE_IO{

public:


    void WriteSimulationParameters(string directory, PARTICLES& Particle, GRID& Grid, SIMULATION_PARAMETERS& SimulationParameters);
    void WriteParticlePositions(string directory, string simulationName, vector<Vector2d>& particlePosition, int frameNumber);
    void WriteParticleVelocities(string directory, string simulationName, vector<Vector2d>& particleVelocity, int frameNumber);
    void WriteCauchyStress(string directory, string simulationName, vector<Matrix2d>& cauchyStress, int frameNumber );
    void WritedPsidF(string directory, string simulationName, vector<Matrix2d>& dpsi_df, int frameNumber);
    void WriteDeformationGradient( string directory, string simulationName, vector<Matrix2d>& deformationGradient, int frameNumber );
    void WriteElasticDeformationGradient( string directory, string simulationName, vector<Matrix2d>& elasticDeformationGradient, int frameNumber );
    void WritePlasticDeformationGradient( string directory,string simulationName,  vector<Matrix2d>& plasticDeformationGradient, int frameNumber );
    void WriteTimePerIteration(string directory, string simulationName, vector<double>& timeStep);
    void WriteTimeStep(string directory, string simulationName, vector<double>& timeStep);


    void WriteOutput(string directory, string simulationName, int frameNumber, vector<Vector2d>& particlePosition, vector<Vector2d>& particleVelocity, vector<Matrix2d>& cauchyStress, vector<Matrix2d>& deformationGradient, vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& plasticDeformationGradient);

    // NEED TO IMPLMENT / DEBUG

    // void ReadParticlePositions(string directory, vector<Vector3d>& particlePosition, int frameNumber);
// 	void ReadParticleVelocities(string directry,vector<Vector3d>& particleVelocity, int frameNumber);
// 	void ReadCauchyStress(string directory, vector<Matrix3d>& cauchyStress, int frameNumber);
// 	void ReadDeformationGradient(string directory, vector<Matrix3d>& deformationGradient, int frameNumber);
// 	void ReadElasticDeformationGradient(string directory, vector<Matrix3d>& elasticDeformationGradient, int frameNumber);
// 	void ReadPlasticDeformationGradient(string directory, vector<Matrix3d>& plasticDeformationGradient, int frameNumber);

};

#endif //MPM_LAB_2D_FILE_IO_H
