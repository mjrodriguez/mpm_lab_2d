//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#ifndef MPM_LAB_2D_GRID_H
#define MPM_LAB_2D_GRID_H


#include <vector>
#include "../Eigen/Dense"

#include "INTERPOLATION.h"
#include "ELASTOPLASTIC.h"

using namespace std;
using namespace Eigen;

class GRID{

private:

    int m_nx;
    int m_ny;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;

    // EQUAL GRID SPACING IN ALL DIRECTIONS
    double m_h;

public:

    /***********************************************************************************************************************
    *
    * THIS IS ALL THE GRID PARAMETERS AND STUFF
    *
    **********************************************************************************************************************/

    vector<double> mass;
    vector<VectorXd> nodeWeights;

    vector<Vector2d> position;
    vector<Vector2d> velocity;
    vector<Vector2d> newVelocity;
    vector<Vector2d> force;

    vector<Vector2i> massList;
    vector<Vector2i> nodeWeightsList;

    void SetDefaultGrid();


    int Index2D(const int i, const int j);

    double GetGridSpacing();

    double Getxmin();

    double Getxmax();

    int GetN();


    void InitializeMassGrid();

    void InitializePositionGrid();

    void InitializeVelocityGrid();

    void InitializeGrid();



    /***********************************************************************************************************************
    *
    *  THIS IS ALL THE GRID COMPUTATION FUNCTIONS
    *
    **********************************************************************************************************************/

    ////////////////////////////////////////////////////////////////////////////////////
    // PARTICLE TO GRID FUNCTION:
    // Inputs: 	massGrid => The mass defined on the grid.
    //			velocityGrid => velocity defined on the grid.
    //			positionParticle => Current position of particles.
    // 			massParticle => mass of the particles.
    //			velocityParticle => Current velocity of particles.
    //
    // Goal of Function: The goal of the function is to transfer the mass and velocity
    //						from the lagrangian state to the Eulerian state.
    ////////////////////////////////////////////////////////////////////////////////////

    void ComputeGridNodeWeights(vector<Vector2d> &positionParticle, INTERPOLATION &Interpolation);


    void
    ParticleToGrid(vector<double> &massParticle, vector<Vector2d> &positionParticle, vector<Vector2d> &velocityParticle,
                   INTERPOLATION &Interpolation);

    void NodesWithMass();

    void
    ComputeGridForces(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<double> &JElastic,
                      vector<double> &JPlastic,
                      vector<double> &volumeParticle, vector<Vector2d> &positionParticle,
                      vector<Matrix2d> &cauchyStress,
                      vector<Matrix2d> &elasticDeformationGradient, vector<Matrix2d> &RE, vector<Matrix2d> &SE,
                      ELASTOPLASTIC &ConstitutiveModel, INTERPOLATION &Interpolation);


    void AddGravityForce(Vector2d &gravity);

    void UpdateVelocityGrid(double timeStep);

    void GridCollision(double frictionCoeff);

    void
    UpdateDeformationGradient(bool usePlasticity, double timeStep, double criticalStretch, double criticalCompression,
                              vector<double> &JElastic,
                              vector<double> &JPlastic, vector<Vector2d> &positionParticle,
                              vector<Matrix2d> &elasticDeformationGradient,
                              vector<Matrix2d> &plasticDeformationGradient, vector<Matrix2d> &deformationGradient,
                              INTERPOLATION &Interpolation);

    void ClearEulerianFields();

    Matrix2d Clamp(Matrix2d &principalValues, double criticalStretch, double criticalCompression);

};

#endif //MPM_LAB_2D_GRID_H