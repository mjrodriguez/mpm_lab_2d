//
// Created by Martin Rodriguez Cruz on 9/19/16.
//

#ifndef MPM_V4_0_PARTICLES_H
#define MPM_V4_0_PARTICLES_H

#include <vector>
#include "../Eigen/Dense"

#include "INTERPOLATION.h"
#include "TOOLS.h"


using namespace std;
using namespace Eigen;


class PARTICLES{

private:
    double m_particlesPerDirection;
    double m_cubeLength;
    double m_numberOfParticles;
    double m_particleGridSpacing;
    Vector2d m_initialVelocity;
    Vector2d m_anchorPoint;
public:
    vector<double> mass;
    vector<double> volume;
    vector<double> density;

    vector<double> JElastic; // Jacobian of Elastic Deformation Gradient
    vector<double> JPlastic; // Jacobian of Plastic Deformaiton Gradient


    //////////////////////////////////////////////////////////////
    // 2D variables
    /////////////////////////////////////////////////////////////

    vector<Vector2d> position;
    vector<Vector2d> velocity;
    vector<Vector2d> velocityPIC;
    vector<Vector2d> velocityFLIP;

    vector< Matrix2d > plasticDeformationGradient;
    vector< Matrix2d > elasticDeformationGradient;
    vector< Matrix2d > deformationGradient;
    vector< Matrix2d > cauchyStress;
    vector< Matrix2d > R; // Rotation Matrix from Polar Decomposition
    vector< Matrix2d > S; // Symmetric Matrix From Polar Decomposition




    /****************************************************
    * SETTING THE PARTICLE STATE
    ****************************************************/
    void SetCube( const double particlesPerDirection, const double cubeLength, const Vector2d anchorPoint );
    void SetDefaultParticles();
    Vector2d GetInitialVelocity();
    void SetInitialVelocity(const Vector2d& v0);
    void InitializeParticles();

    void InitializeParticleMass( );
    void InitializeParticleVolume( );
    void InitializeParticleDensity();
    void InitializeDeformationGradients();
    void InitializeParticlePositions();
    void InitializeParticleVelocities();

    double GetNumberOfParticles();

    /****************************************************
    * UPDATING THE PARTICLE STATE
    ****************************************************/

    void UpdateParticleVelocity( double alpha );
    void UpdateParticlePosition( double timeStep );

    /****************************************************
    * PARTICLE COMPUTATIONS
    ****************************************************/

    void ComputeVolumeDensity(  int N, double h, vector<double>& massGrid, vector<Vector2i>& massList, TOOLS Tools, INTERPOLATION Interpolation );
    void ParticleCollision( double frictionCoeff );
    void ComputePicFlipVelocities( int N, double h, vector<Vector2i>& massList, vector<VectorXd>& gridNodeWeights, vector<Vector2d>& velocityGrid, vector<Vector2d>& newVelocityGrid, TOOLS Tools, INTERPOLATION Interpolation );









};




#endif //MPM_V4_0_PARTICLES_H
