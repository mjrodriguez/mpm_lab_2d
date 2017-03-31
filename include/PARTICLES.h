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

public:
    double m_particlesPerDirection;
    double m_cubeLength;
	double m_area;
    double m_numberOfParticles;
    double m_particleGridSpacing;
	double m_initialDensity;
    Vector2d m_initialVelocity;
    Vector2d m_anchorPoint;
	
	
	string m_particleSimulationName;
	
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
	void SetInitialDensity(const double density);
	
    void SetCube( const double particlesPerDirection, const double cubeLength, const Vector2d anchorPoint );
    void SetRectangle(const double particlesInX, const double particlesInY, const double xLength, const double yLength , const Vector2d anchorPoint);
	// void SetSnowball(const double h, const Vector2d& gridPosition);
	// void Phi(const double x, const double y);
	virtual void SetDefaultParticles() = 0;
    Vector2d GetInitialVelocity();
    void SetInitialVelocity(const Vector2d& v0);
    void InitializeParticles();

    void InitializeParticleMass( );
    void InitializeParticleVolume( );
    void InitializeParticleDensity();
    void InitializeDeformationGradients();
    void InitializeParticlePositions();
    virtual void InitializeParticleVelocities() = 0;

    double GetNumberOfParticles();
	
	string GetParticleSimulationName();

    /****************************************************
    * UPDATING THE PARTICLE STATE
    ****************************************************/

    void UpdateParticleVelocity( double alpha );
    void UpdateParticlePosition( double timeStep );

    /****************************************************
    * PARTICLE COMPUTATIONS
    ****************************************************/
	void ComputeVolumeDensity(  int N, double h, vector<double>& massGrid, vector<Vector2i>& massList, TOOLS Tools, INTERPOLATION Interpolation );
    virtual void ParticleCollision( double frictionCoeff ) = 0;
    void ComputePicFlipVelocities( int N, double h, vector<Vector2i>& massList, vector<VectorXd>& gridNodeWeights, vector<Vector2d>& velocityGrid, vector<Vector2d>& newVelocityGrid, TOOLS Tools, INTERPOLATION Interpolation );









};





/////////////////////////////////////////////////////////////////
//
// CUBE TO CUBE COLLISION
//
//////////////////////////////////////////////////////////////////

class CUBE_TO_CUBE_COLLISION: public PARTICLES {
public:
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff );
};



/////////////////////////////////////////////////////////////////
//
// RECTANGLE FREE FALL with CYLINDER
//
//////////////////////////////////////////////////////////////////

class RECTANGLE_FREEFALL_CYLINDER: public PARTICLES {
public:
	double LevelSet(const double x, const double y);
	Vector2d GradientLevelSet(const double x, const double y);
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff ); 
};

/////////////////////////////////////////////////////////////////
//
// RECTANGLE FREE FALL with HAT OBSTACLE
//
//////////////////////////////////////////////////////////////////

class RECTANGLE_FREEFALL_HAT_OBSTACLE: public PARTICLES {
public:
	double LevelSet(const double x, const double y);
	Vector2d GradientLevelSet(const double x, const double y);
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff ); 
};



/////////////////////////////////////////////////////////////////
//
// CUBE FREE FALL with no OBSTACLE
//
//////////////////////////////////////////////////////////////////

class CUBE_FREEFALL: public PARTICLES {
public:
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff );
};



/////////////////////////////////////////////////////////////////
//
// CUBE CRASH INTO PADDED GROUND
//
//////////////////////////////////////////////////////////////////

class CUBE_CRASH_PADDED_GROUND: public PARTICLES {
public:
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff );
};





/////////////////////////////////////////////////////////////////
//
// CUBE CRASH INTO WALL
//
//////////////////////////////////////////////////////////////////

class CUBE_CRASH_WALL: public PARTICLES {
public:
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff );
};


/////////////////////////////////////////////////////////////////
//
// CUBE_TO_CUBE_FREEFALL
//
//////////////////////////////////////////////////////////////////

class CUBE_TO_CUBE_FREEFALL: public PARTICLES {
public:
	void SetDefaultParticles();
	void InitializeParticleVelocities();
	void ParticleCollision( double frictionCoeff );
};





#endif //MPM_V4_0_PARTICLES_H
