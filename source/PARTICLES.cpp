////////////////////////////////////////////////////////////////
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <vector>
#include "../Eigen/Dense"
#include "../Eigen/Geometry" 
#include "../include/INTERPOLATION.h"
#include "../include/PARTICLES.h"
#include "../include/TOOLS.h"
#include <iostream>
#include <math.h>       /* sin */

#define PI 3.14159265

using namespace std;
using namespace Eigen;



double PARTICLES::GetNumberOfParticles() {
    return m_numberOfParticles;
}

Vector2d PARTICLES::GetInitialVelocity(){
    return m_initialVelocity;
}

string PARTICLES::GetParticleSimulationName(){
	return m_particleSimulationName;
}


////////////////////////////////////////////////////////////////
// SETTING THE PARTICLES STATES
///////////////////////////////////////////////////////////////


void PARTICLES::SetInitialVelocity(const Vector2d &v0) {
    m_initialVelocity = v0;
}


void PARTICLES::SetInitialDensity(const double density){
	m_initialDensity = density;
}




// void PARTICLES::Phi(const double x, const double y){
// 	return (x-0.5)*(x-0.5) + (y - 0.6)*(y-0.6) - 0.1*0.1;
// }
//
// void PARTICLES::SetSnowball(const double h, const Vector2d& gridPosition){
// 	Vector2d P1 = (1 + Vector2d::Random())*0.5*h
// 	Vector2d P2 = (1 + Vector2d::Random())*0.5*h
// 	Vector2d P3 = (1 + Vector2d::Random())*0.5*h
// 	Vector2d P4 = (1 + Vector2d::Random())*0.5*h
//
// 	// FINDING CENTER of CELLS
// 	for (int i = 0; i < gridPosition.size()-1; i++ ){
// 		for (int j = 0; j < gridPosition.size()-1; j++){
// 			double center = gridPosition[]
// 		}
// 	}
//
//
// }




void PARTICLES::InitializeParticles() {

    // ALLOCATING MEMORY
    mass.reserve(m_numberOfParticles);
    volume.reserve(m_numberOfParticles);
    density.reserve(m_numberOfParticles);
    plasticDeformationGradient.reserve(m_numberOfParticles);
    elasticDeformationGradient.reserve(m_numberOfParticles);
    deformationGradient.reserve(m_numberOfParticles);
    cauchyStress.reserve(m_numberOfParticles);
    JElastic.reserve(m_numberOfParticles);
    JPlastic.reserve(m_numberOfParticles);
    R.reserve(m_numberOfParticles);
    S.reserve(m_numberOfParticles);
    velocity.reserve(m_numberOfParticles);
    velocityFLIP.reserve(m_numberOfParticles);
    velocityPIC.reserve(m_numberOfParticles);

    InitializeParticleMass();
    InitializeParticleVolume();
    InitializeParticleDensity();
    InitializeDeformationGradients();
//    InitializeParticlePositions();
    InitializeParticleVelocities();
}



void PARTICLES::SetCube(const double particlesPerDirection, const double cubeLength, const Vector2d anchorPoint){
	
    m_particlesPerDirection = particlesPerDirection;
    m_cubeLength = cubeLength;
	
	m_area += cubeLength*cubeLength;
	
    m_anchorPoint = anchorPoint;
    m_particleGridSpacing = ( m_cubeLength )/m_particlesPerDirection;
    m_numberOfParticles += m_particlesPerDirection*m_particlesPerDirection;
	/*Matrix2d Rotation;
	double THETA = 0;
	Rotation << cos(THETA), -sin(THETA), sin(THETA), cos(THETA);
	cout << "ROTATION MATRIX= " << Rotation << endl;*/
	
    for (int px = 0; px < m_particlesPerDirection; px++) {
        for (int py = 0; py < m_particlesPerDirection; py++) {
			Vector2d P(m_anchorPoint.x() + m_particleGridSpacing*px, m_anchorPoint.y() + m_particleGridSpacing*py );
			// P = Rotation*P;
			// cout  << P << endl;
			
            position.push_back( P );
                // cout << Vector3T (0.25 + m_particleGridSpacing*px, 0.25 + m_particleGridSpacing*py , 0.75 + m_particleGridSpacing*pz) << endl;
        }
    }

}

void PARTICLES::SetRectangle(const double particlesInX, const double particlesInY, const double xLength, const double yLength , const Vector2d anchorPoint){
	
    m_anchorPoint = anchorPoint;
    double xParticleGridSpacing = ( xLength )/particlesInX;
	double yParticleGridSpacing = ( yLength )/particlesInY;
    m_numberOfParticles += particlesInX*particlesInY;
	m_area += xLength*yLength;
	
	
    for (int px = 0; px < particlesInX; px++) {
        for (int py = 0; py < particlesInY; py++) {
            position.push_back( Vector2d (m_anchorPoint.x() + xParticleGridSpacing*px, m_anchorPoint.y() + yParticleGridSpacing*py ) );
                // cout << Vector3T (0.25 + m_particleGridSpacing*px, 0.25 + m_particleGridSpacing*py , 0.75 + m_particleGridSpacing*pz) << endl;
        }
    }

}



void PARTICLES::InitializeParticleMass() {
	
    for (int p = 0; p < m_numberOfParticles; p++){
        mass.push_back((double) (m_initialDensity*m_area)/m_numberOfParticles );
    }
}

void PARTICLES::InitializeParticleVolume() {
    for (int p = 0; p < m_numberOfParticles; p++){
        volume.push_back((double) 0 );
    }
}

void PARTICLES::InitializeParticleDensity() {
    for (int p = 0; p < m_numberOfParticles; p++){
        density.push_back((double) 0 );
    }
}

void PARTICLES::InitializeDeformationGradients() {
    for (int p = 0; p < m_numberOfParticles; p++){
        plasticDeformationGradient.push_back( Matrix2d::Identity() );
        elasticDeformationGradient.push_back( Matrix2d::Identity() );
        deformationGradient.push_back( Matrix2d::Identity() );

        cauchyStress.push_back( Matrix2d::Zero() );

        JElastic.push_back( elasticDeformationGradient[p].determinant() );
        JPlastic.push_back( plasticDeformationGradient[p].determinant() );
        R.push_back( Matrix2d::Identity() );
        S.push_back( Matrix2d::Identity() );

    }
}




////////////////////////////////////////////////////////////////
// UPDATING PARTICLES STATES
///////////////////////////////////////////////////////////




void PARTICLES::UpdateParticleVelocity(double alpha)  {
    for (int p = 0; p < m_numberOfParticles; p++ ){

        velocity[p] = ( (double) 1 - alpha )*velocityPIC[p] + alpha*velocityFLIP[p];

//        velocity[p].x() = ( 1.0 - alpha )*velocityPIC[p].x() + alpha*velocityFLIP[p].x();
//        velocity[p].y() = ( 1.0 - alpha )*velocityPIC[p].y() + alpha*velocityFLIP[p].y();
//        velocity[p].z() = ( 1.0 - alpha )*velocityPIC[p].z() + alpha*velocityFLIP[p].z();
    }
}

void PARTICLES::UpdateParticlePosition(double timeStep) {
    for (int p = 0; p < m_numberOfParticles; p++ ){
        position[p] += timeStep*velocity[p];

//        position[p].x() += timeStep*velocity[p].x();
//        position[p].y() += timeStep*velocity[p].y();
//        position[p].z() += timeStep*velocity[p].z();
    }
}





void PARTICLES::ComputeVolumeDensity( int N, double h, vector<double>& massGrid, vector<Vector2i>& massList, TOOLS Tools, INTERPOLATION Interpolation ) {

    double densityGrid = 0;
    double rho_temp = 0;

    for (int p = 0; p < GetNumberOfParticles(); p++) {
        for (int i = 0; i < massList.size(); i++) {
            densityGrid = massGrid[ Tools.Index2D(N, massList[i].x(), massList[i].y()) ]/(h*h);
            rho_temp += densityGrid*Interpolation.Weight(h, position[p], massList[i].x(), massList[i].y()  );
        }
        density[p] = rho_temp;
        volume[p] = mass[p]/rho_temp;

        densityGrid = 0;
        rho_temp = 0;

    }

}



void PARTICLES::ComputePicFlipVelocities(int N, double h, vector<Vector2i>& massList,  vector<VectorXd>& gridNodeWeights, vector<Vector2d>& velocityGrid, vector<Vector2d>& newVelocityGrid, TOOLS Tools, INTERPOLATION Interpolation ) {
    for (int p = 0; p < m_numberOfParticles; p++) {
        velocityFLIP[p] = velocity[p]; // Vector3d::Zero();
        velocityPIC[p] = Vector2d::Zero();



        for (int i = 0; i < massList.size(); i++) {
            int ig, jg;
            ig = massList[i].x();
            jg = massList[i].y();
            double wip = gridNodeWeights[i][p]; //Interpolation.Weight(h, position[p], ig,jg,kg);
            velocityPIC[p] += newVelocityGrid[ Tools.Index2D(N,ig,jg) ]*wip;
            velocityFLIP[p] += ( newVelocityGrid[Tools.Index2D(N,ig,jg)] - velocityGrid[ Tools.Index2D(N,ig,jg) ] )*wip;
        }

        // cout << "Particle " << p << ": Difference of PIC and FLIP =  " << (velocityPIC[p].transpose() - velocityFLIP[p].transpose()).norm() << endl;

        // velocityFLIP[p] += velocity[p];

    }

}










/////////////////////////////////////////////////////////////////
//
// CUBE TO CUBE COLLISION:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


void CUBE_TO_CUBE_COLLISION::SetDefaultParticles() {
	
	m_particleSimulationName = string("cube_to_cube_collision");
	
	m_numberOfParticles = 0;
	m_area = 0;
	
    Vector2d anchor1(0.01,0.5);
    SetCube( 25, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	Vector2d anchor2(0.79, 0.4);
    SetCube( 25, 0.2, anchor2 );

    InitializeParticles();
}


void CUBE_TO_CUBE_COLLISION::InitializeParticleVelocities() {
	Vector2d m_initialVelocity = Vector2d(5,3);
    for (int p = 0; p < m_numberOfParticles/2; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }
	
	for (int p = m_numberOfParticles/2; p < m_numberOfParticles; p++){
		// Collision Against Each Other
        velocity.push_back(Vector2d(-m_initialVelocity.x(), m_initialVelocity.y()));
    	velocityPIC.push_back( Vector2d( 0, 0) );
    	velocityFLIP.push_back( Vector2d( 0, 0) );
	}
}




void CUBE_TO_CUBE_COLLISION::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;

		
		///////////////////////////////////////
		// GROUND COLLISION
		///////////////////////////////////////
	
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
        }
        else if (position[p].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
        }
        else if (position[p].x() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}





/////////////////////////////////////////////////////////////////
//
// RECTANGLE_FREEFALL_OBSTACLE:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


double RECTANGLE_FREEFALL_CYLINDER::LevelSet(const double x, const double y){
	return (x-0.5)*(x-0.5) + (y-0.3)*(y-0.3) - 0.1*0.1;
}
Vector2d RECTANGLE_FREEFALL_CYLINDER::GradientLevelSet(const double x, const double y){
	Vector2d v = Vector2d( 2*(x-0.5), 2*(y-0.3) );
	return v/v.norm();
}


void RECTANGLE_FREEFALL_CYLINDER::SetDefaultParticles() {
	
	m_particleSimulationName = string("rectangle_freefall_cylinder");
	
	m_numberOfParticles = 0;
	m_area = 0;
	
	/*
    Vector2d anchor1(0.1,0.5);
    SetCube( 20, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	Vector2d anchor2(0.7, 0.5)
    SetCube( 20, 0.2, anchor2 );
	*/
	
    Vector2d anchor1(0.3,0.75);
    SetRectangle( 50, 25, 0.4, 0.2, anchor1 );
	

    InitializeParticles();
}


void RECTANGLE_FREEFALL_CYLINDER::InitializeParticleVelocities() {
	m_initialVelocity = Vector2d(0,-3);
    for (int p = 0; p < m_numberOfParticles; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
	}
}


void RECTANGLE_FREEFALL_CYLINDER::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;


		//////////////////////////////////////////
		// COLLISION WITH CYLINDER
		//////////////////////////////////////////

		double phi = LevelSet(position[p].x(),position[p].y());
		if ( phi <= (double) 0 ){
			hasItCollided = true;
			collisionNormal = GradientLevelSet( position[p].x(), position[p].y() );
		}
		
		// GROUND COLLISION
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
        }
        else if (position[p].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
        }
        else if (position[p].x() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}


/////////////////////////////////////////////////////////////////
//
// RECTANGLE_FREEFALL_HAT_OBSTACLE:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


double RECTANGLE_FREEFALL_HAT_OBSTACLE::LevelSet(const double x, const double y){
	if ( x >= 0.4 && x <= 0.6 && y >= 0.5 && y <= 0.6 ){
		return -( 0.6 - fabs(x-0.5)-y );
	}	
	else{
		return 1;
	}
}
Vector2d RECTANGLE_FREEFALL_HAT_OBSTACLE::GradientLevelSet(const double x, const double y){
	Vector2d v;
	if ( x >= 0.4 && x <= 0.5 && y >= 0.5 && y <= 0.6 ){
		v = Vector2d( -1, 1 );
	}
	else if ( x > 0.5 && x <= 0.6 && y >= 0.5 && y <= 0.6 ){
		v = Vector2d( 1, 1 );
	}
	return v/v.norm();
}


void RECTANGLE_FREEFALL_HAT_OBSTACLE::SetDefaultParticles() {
	m_particleSimulationName = string("rectangle_freefall_hat_obstacle");
	m_numberOfParticles = 0;
	m_area = 0;
	
	/*
    Vector2d anchor1(0.1,0.5);
    SetCube( 20, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	Vector2d anchor2(0.7, 0.5)
    SetCube( 20, 0.2, anchor2 );
	*/
	
    Vector2d anchor1(0.3,0.75);
    SetRectangle( 50, 25, 0.4, 0.2, anchor1 );
	

    InitializeParticles();
}


void RECTANGLE_FREEFALL_HAT_OBSTACLE::InitializeParticleVelocities() {
	m_initialVelocity = Vector2d(0,-3);
    for (int p = 0; p < m_numberOfParticles; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
	}
}


void RECTANGLE_FREEFALL_HAT_OBSTACLE::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0, 0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;


		//////////////////////////////////////////
		// COLLISION WITH CYLINDER
		//////////////////////////////////////////

		double phi = LevelSet(position[p].x(),position[p].y());

		if ( phi <= (double) 0 ){
			// cout << "Collision with Cylinder PARTICLE" << endl;
			hasItCollided = true;
			collisionNormal = GradientLevelSet( position[p].x(), position[p].y() );
		}
		// GROUND COLLISION
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
        }
        else if (position[p].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
        }
        else if (position[p].x() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}



/////////////////////////////////////////////////////////////////
//
// CUBE TO CUBE COLLISION:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


void CUBE_CRASH_WALL::SetDefaultParticles() {
	stickyCollision = true;
	m_particleSimulationName = string("cube_crash_wall");
	
	m_numberOfParticles = 0;
	m_area = 0;
	
    Vector2d anchor1(0.2,0.4);
    SetCube( 40, 0.15, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	// Vector2d anchor2(0.7, 0.5)
    // SetCube( 20, 0.2, anchor2 );

    InitializeParticles();
}


void CUBE_CRASH_WALL::InitializeParticleVelocities() {
	Vector2d m_initialVelocity = Vector2d(3,3.0);
    for (int p = 0; p < m_numberOfParticles; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d( m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }
}




void CUBE_CRASH_WALL::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;

		
		///////////////////////////////////////
		// GROUND COLLISION
		///////////////////////////////////////
	
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
			stickyCollision = true;
        }
        else if (position[p].x() <= (double) 0.05 ){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
			stickyCollision = true;
        }
        else if (position[p].x() >= (double) 0.95 ){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
			stickyCollision = true;
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity || stickyCollision){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}









/////////////////////////////////////////////////////////////////
//
// CUBE_CRASH_PADDED_GROUND:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


void CUBE_CRASH_PADDED_GROUND::SetDefaultParticles() {
	
	m_particleSimulationName = string("cube_crash_padded_ground");
	
	m_numberOfParticles = 0;
	m_area = 0;
	
    Vector2d anchor1(0.1,0.5);
	double cubeParticlesSide = 30;
    SetCube( cubeParticlesSide, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	// Vector2d anchor2(0.7, 0.5)
    // SetCube( 20, 0.2, anchor2 );
	
	
    Vector2d anchor2(0.05,0.05);
    SetRectangle( 250, 20, 0.9, 0.1, anchor2 );
	
    InitializeParticles();
}


void CUBE_CRASH_PADDED_GROUND::InitializeParticleVelocities() {
	Vector2d m_initialVelocity = Vector2d(2,-2);
    for (int p = 0; p < 30*30; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }
    for (int p = 30*30; p < m_numberOfParticles; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(0, 0));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }

}

void CUBE_CRASH_PADDED_GROUND::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;

		
		///////////////////////////////////////
		// GROUND COLLISION
		///////////////////////////////////////
	
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
        }
        else if (position[p].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
        }
        else if (position[p].x() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}










/////////////////////////////////////////////////////////////////
//
// CUBE TO CUBE COLLISION:
// INHERIT THE MAIN FUNCTIONS FROM PARTICLES CLASS
// SETS UP DIFFERENT SIMULATION CASES
//
//////////////////////////////////////////////////////////////////


void CUBE_TO_CUBE_FREEFALL::SetDefaultParticles() {
	
	m_particleSimulationName = string("cube_to_cube_freefall");
	
	m_numberOfParticles = 0;
	m_area = 0;
	
    Vector2d anchor1(0.4,0.5);
    SetCube( 30, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	Vector2d anchor2(0.5, 0.05);
	SetCube( 30, 0.2, anchor2 );

    InitializeParticles();
}


void CUBE_TO_CUBE_FREEFALL::InitializeParticleVelocities() {
	Vector2d m_initialVelocity = Vector2d(0,-2);
    for (int p = 0; p < m_numberOfParticles/2; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }
    for (int p = m_numberOfParticles/2; p < m_numberOfParticles; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(0, 0));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }

}




void CUBE_TO_CUBE_FREEFALL::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;
		
		///////////////////////////////////////
		// GROUND COLLISION
		///////////////////////////////////////
	
        if (position[p].y() <= ( (double) 0.05 ) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,1);
        }
        else if (position[p].y() >= ( (double) 0.95) ){
            hasItCollided = true;
            collisionNormal = Vector2d(0,-1);
        }
        else if (position[p].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1,0);
        }
        else if (position[p].x() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(-1,0);
        }

        if (hasItCollided){
            Vector2d velocityRel;
            velocityRel = velocity[p] - velocityCollision;
            double normalVelocity = velocityRel.dot(collisionNormal);


            if (normalVelocity < 0){
                Vector2d tangentVelocity;
                tangentVelocity = velocityRel - collisionNormal * normalVelocity;


                Vector2d velocityRelPrime;
                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity){
                    velocityRelPrime = Vector2d(0,0);
                }
                else {
                    double one_over_vt_norm = 1/tangentVelocity.norm();
                    velocityRelPrime = tangentVelocity + frictionCoeff*normalVelocity*one_over_vt_norm*tangentVelocity;
                }

                velocity[p] = velocityRelPrime + velocityCollision;

            }



        }
    }

}
