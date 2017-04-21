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

public:

    int m_nx;
    int m_ny;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;

    // EQUAL GRID SPACING IN ALL DIRECTIONS
    double m_h;
	
	

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
	double Getymin();
    double Getxmax();

    int GetN();


    void InitializeMassGrid();

    void InitializePositionGrid();

    void InitializeVelocityGrid();
	void InitializeForceGrid();
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
	void ComputeNeighborhood(Vector2d& positionParticle, Vector2i& gridIndexCenter, Vector2i& minGridIndex, Vector2i& maxGridIndex);

    void ParticleToGrid(vector<double> &massParticle, vector<Vector2d> &positionParticle, vector<Vector2d> &velocityParticle,
                   INTERPOLATION &Interpolation);

    void NodesWithMass();

    void ComputeGridForces(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<double> &JElastic,
                      vector<double> &JPlastic,
                      vector<double> &volumeParticle, vector<Vector2d> &positionParticle,
                      vector<Matrix2d> &cauchyStress,
                      vector<Matrix2d> &elasticDeformationGradient, vector<Matrix2d> &RE, vector<Matrix2d> &SE,
                      ELASTOPLASTIC &ConstitutiveModel, INTERPOLATION &Interpolation);


    void AddGravityForce(Vector2d &gravity);

    void UpdateVelocityGrid(double timeStep);
	
	void UpdateGravityVelocityGrid(const double timeStep, const Vector2d& gravity);

	////////////////////////////////////////////////////////
	// LEVEL SETS and COLLISIONS
	/////////////////////////////////////////////////////////
	virtual void GridCollision(double frictionCoeff) = 0;

    void UpdateDeformationGradient(bool usePlasticity, double timeStep, double criticalStretch, double criticalCompression,
                              vector<double> &JElastic,
                              vector<double> &JPlastic, vector<Vector2d> &positionParticle,
                              vector<Matrix2d> &elasticDeformationGradient,
                              vector<Matrix2d> &plasticDeformationGradient, vector<Matrix2d> &deformationGradient,
                              INTERPOLATION &Interpolation);

    void ClearEulerianFields();

    Matrix2d Clamp(Matrix2d &principalValues, double criticalStretch, double criticalCompression);



    /***********************************************************************************************************************
    *
    *  IMPLICIT TIME INTEGRATION MODULES
    *
    **********************************************************************************************************************/
	
	
	void ImplicitUpdateVelocityGrid(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);
	void ConjugateResidual(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);
	vector<Vector2d> ComputeHessianAction(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);
	vector<Matrix2d> ComputeAp(const double timeStep, const bool usePlasticity, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<Vector2d>& positionParticle, vector<Vector2d>& u, const vector<Matrix2d>& F, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);
	vector<Vector2d> ComputeDeltaForce(const vector<double>& particleVolume, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& F, const vector<Matrix2d>& Ap, INTERPOLATION& Interpolation );
	// vector<Matrix2d> ComputeDeltaF(const double timeStep, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& elasticDeformationGradient, INTERPOLATION& Interpolation);
	// void ComputeFHat(double dt, vector<double>& JE_hat, const vector<Vector2d>& particlePosition, const vector<Matrix2d>& F, vector<Matrix2d>& F_hat, vector<Matrix2d>&  R_hat, vector<Matrix2d>&  S_hat, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);
	void ComputeFHat( const double dt, const Matrix2d& velocityGradient, vector<double>& JE_hat, const vector<Vector2d>& particlePosition, const Matrix2d& F, vector<Matrix2d>& F_hat, vector<Matrix2d>&  R_hat, vector<Matrix2d>&  S_hat, INTERPOLATION& Interpolation);
	vector<Matrix2d> ComputeDeltaF(const double timeStep, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, vector<double>& JE_hat, const vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& F_hat, vector<Matrix2d>& R_hat, vector<Matrix2d>& S_hat , INTERPOLATION& Interpolation);
	vector<Matrix2d> ComputeDeltaJFiT(const vector<Matrix2d>& F, const vector<Matrix2d>& deltaF);
	vector<Matrix2d> Compute_dJFiT(const Matrix2d& Fe);
	vector<Matrix2d> ComputeDeltaR(const vector<Matrix2d>& deltaF, const vector<Matrix2d>& S, const vector<Matrix2d>& R);
	Matrix2d DoubleDotProduct(const vector<Matrix2d>& C, const Matrix2d& B);
	double DoubleDotProduct(const Matrix2d& C, const Matrix2d& B);
	double InnerProduct(vector<Vector2d>& a, vector<Vector2d>& b);
	vector<Vector2d> scalarVectorProduct(double a, vector<Vector2d>& p);
	vector<Vector2d> VectorAddition(vector<Vector2d>& p, vector<Vector2d>& q);
	vector<Vector2d> VectorSubtraction(vector<Vector2d>& p, vector<Vector2d>& q);
	
};


class GRID_NO_OBSTACLE: public GRID {
public:
	void GridCollision(double frictionCoeff);
};

class GRID_CYLINDER_OBSTACLE: public GRID{
public:
	double LevelSet(const double x, const double y);
	Vector2d GradientLevelSet(const double x, const double y);
	void GridCollision(double frictionCoeff);
};


class GRID_HAT_OBSTACLE: public GRID{
public:
	double LevelSet(const double x, const double y);
	Vector2d GradientLevelSet(const double x, const double y);
	void GridCollision(double frictionCoeff);
};


#endif //MPM_LAB_2D_GRID_H