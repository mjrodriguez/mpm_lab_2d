////////////////////////////////////////////////////////////////
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <vector>
#include "../Eigen/Dense"

#include "../include/INTERPOLATION.h"
#include "../include/PARTICLES.h"
#include "../include/TOOLS.h"


using namespace std;
using namespace Eigen;



double PARTICLES::GetNumberOfParticles() {
    return m_numberOfParticles;
}

Vector2d PARTICLES::GetInitialVelocity(){
    return m_initialVelocity;
}



////////////////////////////////////////////////////////////////
// SETTING THE PARTICLES STATES
///////////////////////////////////////////////////////////////




void PARTICLES::SetCube(const double particlesPerDirection, const double cubeLength, const Vector2d anchorPoint){
	
    m_particlesPerDirection = particlesPerDirection;
    m_cubeLength = cubeLength;
    m_anchorPoint = anchorPoint;
    m_particleGridSpacing = ( m_cubeLength )/m_particlesPerDirection;
    m_numberOfParticles += m_particlesPerDirection*m_particlesPerDirection;

    for (int px = 0; px < m_particlesPerDirection; px++) {
        for (int py = 0; py < m_particlesPerDirection; py++) {
            position.push_back( Vector2d (m_anchorPoint.x() + m_particleGridSpacing*px, m_anchorPoint.y() + m_particleGridSpacing*py ) );
                // cout << Vector3T (0.25 + m_particleGridSpacing*px, 0.25 + m_particleGridSpacing*py , 0.75 + m_particleGridSpacing*pz) << endl;
        }
    }

}



void PARTICLES::SetDefaultParticles() {
	m_numberOfParticles = 0;
	
    Vector2d anchor1(0.3,0.5);
    SetCube( 20, 0.2, anchor1 );
    // Vector2d anchor2(0.7, 0.5);
	Vector2d anchor2(0.2, 0.0)
    SetCube( 20, 0.2, anchor2 );

    InitializeParticles();
}


void PARTICLES::SetInitialVelocity(const Vector2d &v0) {
    m_initialVelocity = v0;
}



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



void PARTICLES::InitializeParticleMass() {
    for (int p = 0; p < m_numberOfParticles; p++){
        mass.push_back((double) (400*0.4)/m_numberOfParticles );
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

void PARTICLES::InitializeParticleVelocities() {
    for (int p = 0; p < m_numberOfParticles/2; p++ ){
        //velocity.push_back( Vector2d( 0.0, 5.0 ) );
        velocity.push_back(Vector2d(m_initialVelocity.x(), m_initialVelocity.y()));
        velocityPIC.push_back( Vector2d( 0, 0) );
        velocityFLIP.push_back( Vector2d( 0, 0) );
    }
	
	for (int p = m_numberOfParticles/2; p < m_numberOfParticles; p++){
		// Free Fall
		velocity.push_back( Vector2d(0, 0) );
		
		// Collision Against Each Other
        // velocity.push_back(Vector2d(-m_initialVelocity.x(), m_initialVelocity.y()));
    	velocityPIC.push_back( Vector2d( 0, 0) );
    	velocityFLIP.push_back( Vector2d( 0, 0) );
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



void PARTICLES::ParticleCollision(double frictionCoeff) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d(0,0);


    for (int p = 0; p < m_numberOfParticles; p++){
        bool hasItCollided = false;

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










