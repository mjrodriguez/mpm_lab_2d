//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <vector>
#include "../Eigen/Dense"

#include "../include/GRID.h"
#include "../include/INTERPOLATION.h"
#include "../include/ELASTOPLASTIC.h"
#include <iostream>

using namespace std;
using namespace Eigen;


/***********************************************************************************************************************
*
* THIS IS ALL THE GRID PARAMETERS AND STUFF
*
**********************************************************************************************************************/


void GRID::SetDefaultGrid() {
    m_nx = 64;
    m_ny = m_nx;

    m_xmin = (double) 0;
    m_xmax = (double) 1;
    m_ymin = m_xmin;
    m_ymax = m_xmax;


    m_h = (m_xmax - m_xmin)/(m_nx - 1);


    mass.reserve(m_nx*m_ny);
    position.reserve(m_nx*m_ny);
    velocity.reserve(m_nx*m_ny);
    newVelocity.reserve(m_nx*m_ny);

    InitializeGrid();
}

int GRID::Index2D(const int i, const int j) {
    return i*m_ny + j;
}

double GRID::GetGridSpacing() {
    return m_h;
}

double GRID::Getxmin() {
    return m_xmin;
}

double GRID::Getxmax() {
    return m_xmax;
}

int GRID::GetN() {
    return m_nx;
}


void GRID::InitializeMassGrid() {
    for (int i = 0; i < m_nx*m_ny; i++){
        mass.push_back((double) 0);
    }
}

void GRID::InitializePositionGrid() {
    for (int i = 0; i < m_nx; i++){
        for (int j = 0; j < m_ny; j++){
            position.push_back( Vector2d( i*m_h, j*m_h) );
        }
    }
}

void GRID::InitializeVelocityGrid() {
    for (int i = 0; i < m_nx*m_ny; i++){
        velocity.push_back( Vector2d (0, 0) );
        newVelocity.push_back( Vector2d (0, 0) );
    }
}


void GRID::InitializeGrid() {
    InitializePositionGrid();
    InitializeVelocityGrid();
    InitializeMassGrid();
}


/***********************************************************************************************************************
*
*  THIS IS ALL THE GRID COMPUTATION FUNCTIONS
*
**********************************************************************************************************************/

void GRID::ComputeGridNodeWeights( vector<Vector2d>& positionParticle, INTERPOLATION& Interpolation ){
    VectorXd wip;
    wip = VectorXd::Zero(positionParticle.size(),1);

    for (int i = 0; i < GetN(); i++){
        for (int j = 0; j < GetN(); j++){

                for (int p = 0; p < positionParticle.size(); p ++){
                    wip[p] = Interpolation.Weight(m_h, positionParticle[p], i, j);
                }

                if (wip.norm() != 0){
                    nodeWeights.push_back(wip);
                    nodeWeightsList.push_back( Vector2i(i,j) );
                }

                wip = VectorXd::Zero(positionParticle.size(),1);



        }
    }
	
	
	
	/*for (int p = 0; p < positionParticle.size(); p++){
		
		for (int i = 0; i < GetN(); i++){
			for (int j = 0; j < GetN(); j++){
				wip[p] = Interpolation.Weight(m_h, positionParticle[p], i, j);
			}
		}
		
		
		
		
	}*/
	
}



void GRID::ParticleToGrid( vector<double>& massParticle, vector<Vector2d>& positionParticle, vector<Vector2d>& velocityParticle, INTERPOLATION& Interpolation) {

    double massWeighted = 0;
    Vector2d momentumWeighted = Vector2d (0, 0);


    for (int i = 0; i < nodeWeightsList.size(); i++){

        // Setting the node indices to something easier to type
        int iw = nodeWeightsList[i].x();
        int jw = nodeWeightsList[i].y();

        for (int p = 0; p < massParticle.size(); p++){

            massWeighted += massParticle[p]*nodeWeights[i][p];
            momentumWeighted.x() += velocityParticle[p].x()*massParticle[p]*nodeWeights[i][p];
            momentumWeighted.y() += velocityParticle[p].y()*massParticle[p]*nodeWeights[i][p];
        }


        // cout << "Weighted Mass = " << massWeighted << endl;

        mass[ Index2D(iw, jw) ] = massWeighted;
        if (mass[ Index2D(iw, jw)] != 0 ){
            velocity[Index2D(iw, jw)].x() = momentumWeighted.x()/mass[Index2D(iw, jw)];
            velocity[Index2D(iw, jw)].y() = momentumWeighted.y()/mass[Index2D(iw, jw)];
        }

        massWeighted = 0;
        momentumWeighted = Vector2d (0, 0);



    }


    /* double wip = 0;
     double massWeighted = 0;
     Vector3d momentumWeighted = Vector3d (0, 0, 0);


     for (int i = 0; i < GetN(); i++){
         for (int j = 0; j < GetN(); j++){
             for (int k = 0; k < GetN(); k++){

                 for (int p = 0; p < massParticle.size(); p++){

                     wip = Interpolation.Weight(m_h, positionParticle[p], i, j, k);
                     massWeighted += massParticle[p]*wip;
                     momentumWeighted.x() += velocityParticle[p].x()*massParticle[p]*wip;
                     momentumWeighted.y() += velocityParticle[p].y()*massParticle[p]*wip;
                     momentumWeighted.z() += velocityParticle[p].z()*massParticle[p]*wip;


                 }


                 mass[Index3D(i,j,k)] = massWeighted;
                 if (mass[Index3D(i,j,k)] != 0 ){
                     velocity[Index3D(i,j,k)].x() = momentumWeighted.x()/mass[Index3D(i,j,k)];
                     velocity[Index3D(i,j,k)].y() = momentumWeighted.y()/mass[Index3D(i,j,k)];
                     velocity[Index3D(i,j,k)].z() = momentumWeighted.z()/mass[Index3D(i,j,k)];
                 }

                 wip = 0;
                 massWeighted = 0;
                 momentumWeighted = Vector3d (0, 0, 0);


             }
         }
     }*/








}


void GRID::NodesWithMass() {
    /*for (int i = 0; i < GetN(); i++){
        for (int j = 0; j < GetN(); j++){
            for (int k = 0; k < GetN(); k++){

                if (mass[Index3D(i,j,k)] != 0 ){
                    massList.push_back(Vector3i(i,j,k));
                }


            }
        }
    }*/

    massList = nodeWeightsList;

    cout << "Length of massList = " << massList.size() << endl;
    cout << "Length of nodeWeightsList = " << nodeWeightsList.size() << endl;

    // cin.get();

}


void GRID::ComputeGridForces(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<double>& JElastic, vector<double>& JPlastic,
                             vector<double>& volumeParticle, vector<Vector2d>& positionParticle, vector<Matrix2d>& cauchyStress,
                             vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& RE, vector<Matrix2d>& SE, ELASTOPLASTIC& ConstitutiveModel, INTERPOLATION& Interpolation) {

    force.reserve(massList.size());

    ConstitutiveModel.CauchyStress(usePlasticity, mu0, lambda0, hardeningCoeff, cauchyStress, JElastic, JPlastic, elasticDeformationGradient, RE, SE);

//    cout << "Cauchy Stress Value = " << Tools.MaxNormValue(Particle.cauchyStress) << endl;
    // cout << SIGMA[1] << endl;

    for (int i = 0; i < massList.size(); i++){
        int ig, jg;

        ig = massList[i].x();
        jg = massList[i].y();


        Vector2d FORCE;
        FORCE = Vector2d::Zero();

        for (int p = 0; p < JElastic.size(); p++){
            FORCE += -JElastic[p]*JPlastic[p]*volumeParticle[p]*cauchyStress[p]*Interpolation.GradientWeight(m_h, positionParticle[p], ig, jg);
        }


        force.push_back(FORCE);

    } // END MASS LIST FOR LOOP
}



void GRID::AddGravityForce(Vector2d& gravity) {
    for (int i = 0; i < massList.size(); i++){
        int ig, jg;

        ig = massList[i].x();
        jg = massList[i].y();


        force[i].x() += mass[ Index2D(ig, jg) ]*gravity.x();
        force[i].y() += mass[ Index2D(ig, jg) ]*gravity.y();

    }
}


void GRID::UpdateVelocityGrid(double timeStep) {
    for (int i = 0; i < massList.size(); i++){
        int ig, jg;

        ig = massList[i].x();
        jg = massList[i].y();

        newVelocity[ Index2D(ig, jg) ].x() = velocity[ Index2D(ig, jg) ].x() + timeStep*force[i].x()/mass[ Index2D(ig, jg) ];
        newVelocity[ Index2D(ig, jg) ].y() = velocity[ Index2D(ig, jg) ].y() + timeStep*force[i].y()/mass[ Index2D(ig, jg) ];

    }
}






///////////////////////////////////////////////////////////////////
//
// LEVEL SETS AND COLLISIONS
//
////////////////////////////////////////////////////////////////////







void GRID::UpdateDeformationGradient(bool usePlasticity, double timeStep, double criticalStretch, double criticalCompression, vector<double>& JElastic,
                                     vector<double>& JPlastic, vector<Vector2d>& positionParticle, vector<Matrix2d>& elasticDeformationGradient,
                                     vector<Matrix2d>& plasticDeformationGradient, vector<Matrix2d>& deformationGradient,
                                     INTERPOLATION& Interpolation) {


    double dt = timeStep;

    for (int p = 0; p < positionParticle.size(); p++){

        Vector2d gradientWeight;
        Matrix2d velocity_Gradient, A;
        velocity_Gradient = Matrix2d::Zero();
        A = Matrix2d::Zero();

        for ( int i = 0; i < massList.size(); i++ ) {
            int ig, jg;
            ig = massList[i].x();
            jg = massList[i].y();
            gradientWeight = Interpolation.GradientWeight(m_h, positionParticle[p], ig, jg);
            A = velocity[ Index2D(ig, jg) ]*gradientWeight.transpose();


            //                cout << "Iteration = " << i << endl;
            //                cout << " Velocity Vector: " << endl;
            //                cout << velocityGrid[ GridNodes.Index3D(ig, jg, kg) ] << endl;
            //
            //                cout << "Gradient Weight: " << endl;
            //                cout << gradientWeight.transpose() << endl;

            velocity_Gradient += A  ;

            //                cout << "Check :\n " << velocityGrid[ GridNodes.Index3D(ig, jg, kg) ]*gradientWeight.transpose() << endl;
            //                cout << "Gradient Velocity : " << endl;
            //                cout << velocity_Gradient << endl;
            //                cout << " " << endl;
        } // END WEIGHT INTERPOLATION TRANSFER

        Matrix2d FE_hat, FP_hat, F_hat;
        FE_hat = Matrix2d::Zero();
        FP_hat = Matrix2d::Zero();
        F_hat = Matrix2d::Zero();

        FE_hat = elasticDeformationGradient[p] + dt*velocity_Gradient*elasticDeformationGradient[p];
        FP_hat = plasticDeformationGradient[p];

        F_hat = FE_hat*FP_hat;

        JacobiSVD<Matrix2d> svd(FE_hat, ComputeFullU | ComputeFullV);
        Matrix2d sigmaMatrix;


        sigmaMatrix = svd.singularValues().asDiagonal();
        if ( usePlasticity == true ){
            sigmaMatrix = Clamp( sigmaMatrix, criticalStretch, criticalCompression );
            elasticDeformationGradient[p] = svd.matrixU()*sigmaMatrix*svd.matrixV().transpose();
            plasticDeformationGradient[p] = svd.matrixV()*sigmaMatrix.inverse()*svd.matrixU().transpose()*F_hat;
        }
        else{
            elasticDeformationGradient[p] = svd.matrixU()*sigmaMatrix*svd.matrixV().transpose();
            plasticDeformationGradient[p] = Matrix2d::Identity(); // svd.matrixV()*sigmaMatrix.inverse()*svd.matrixU().transpose()*F_hat;
        }


        /*cout << "Elastic Deformaiton Gradient = \n" <<  elasticDeformationGradient[p] << endl;
        cout << "Platic Deformation Gradient = \n" << plasticDeformationGradient[p] << endl;
        cout << "Deformation Gradient = \n" << elasticDeformationGradient[p]*plasticDeformationGradient[p] << endl;*/

        deformationGradient[p] = elasticDeformationGradient[p]*plasticDeformationGradient[p];
        JElastic[p] = elasticDeformationGradient[p].determinant();
        JPlastic[p] = plasticDeformationGradient[p].determinant();


    } // END PARTICLE LOOP


}




void GRID::ClearEulerianFields() {
    massList.erase( massList.begin(), massList.end() );
    force.erase( force.begin(), force.end() );
    nodeWeights.erase(nodeWeights.begin(), nodeWeights.end() );
    nodeWeightsList.erase(nodeWeightsList.begin(),nodeWeightsList.end() );

    for (int i = 0; i < velocity.size(); i ++ ){
        mass[i] = (double) 0;
        velocity[i] = Vector2d (0, 0 );
        newVelocity[i] = Vector2d (0, 0);
    }
}

Matrix2d GRID::Clamp(Matrix2d &principalValues, double criticalStretch, double criticalCompression) {
    for (int i = 0; i < 2; i++){

        if ( principalValues(i,i) < 1-criticalCompression ){
            principalValues(i,i) = 1-criticalCompression;
        }
        else if ( principalValues(i,i) > 1+criticalStretch ){
            principalValues(i,i) = 1+criticalStretch;
        }



    }
    Matrix2d A;
    A = principalValues;
    return A;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////
//====================================================================================================//
// GRID_NO_OBSTACLE
//====================================================================================================//
/////////////////////////////////////////////////////////////////////////////////////////////////////////






void GRID_NO_OBSTACLE::GridCollision( double frictionCoeff ) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d ( 0, 0);

    for (int i = 0; i < massList.size(); i++){

        bool hasItCollided = false;

        int ig, jg;
        ig = massList[i].x();
        jg = massList[i].y();
		
		

        // GROUND COLLISION
        if ( position[Index2D(ig,jg)].y() <= (double) 0.05 ){
            hasItCollided = true;
            collisionNormal = Vector2d( 0, 1);
        }
        else if(position[Index2D(ig,jg)].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1, 0);
        }
        else if( position[Index2D(ig,jg)].y() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(0, -1.0);
        }
        else if( position [Index2D(ig,jg)].x() >= (double) 0.95 ){
            hasItCollided = true;
            collisionNormal = Vector2d(-1.0, 0);
        }



        if ( hasItCollided ){

            Vector2d velocityRel = velocity[Index2D(ig,jg)] - velocityCollision;

            double normalVelocity = velocityRel.dot(collisionNormal);

            if (normalVelocity < 0){

                Vector2d tangentVelocity = velocityRel - collisionNormal*normalVelocity;
                Vector2d velocityRelPrime;

                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity ){
                    velocityRelPrime = Vector2d (0, 0);
                }
                else{
                    velocityRelPrime.x() = tangentVelocity.x() + frictionCoeff*normalVelocity*tangentVelocity.x()/tangentVelocity.norm();
                    velocityRelPrime.y() = tangentVelocity.y() + frictionCoeff*normalVelocity*tangentVelocity.y()/tangentVelocity.norm();
            }

                velocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

            } // if normalVelocity < 0

        } // if hasItCollided




        //            cout << "Grid Node ( " << ig << " , " <<  jg << " , " << kg << " ) = " <<  hasItCollided << endl;

    } // for massList

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//====================================================================================================//
// GRID_CYLINDER_OBSTACLE
//====================================================================================================//
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/// CYLINDER
double GRID_CYLINDER_OBSTACLE::LevelSet(const double x, const double y){
	return (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.1*0.1;
}
Vector2d GRID_CYLINDER_OBSTACLE::GradientLevelSet(const double x, const double y){
	Vector2d v = Vector2d( 2*(x-0.5), 2*(y-0.5) );
	return v/v.norm();
}

void GRID_CYLINDER_OBSTACLE::GridCollision( double frictionCoeff ) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d ( 0, 0);

    for (int i = 0; i < massList.size(); i++){

        bool hasItCollided = false;

        int ig, jg;
        ig = massList[i].x();
        jg = massList[i].y();
		
		
		
		//////////////////////////////////////////
		// COLLISION WITH CYLINDER
		//////////////////////////////////////////
		double phi = LevelSet( position[Index2D(ig,jg)].x(), position[Index2D(ig,jg)].y() );
		if ( phi <= 0 ){
			// cout << "Collision with Cylinder GRID" << endl;
			hasItCollided = true;
			collisionNormal = GradientLevelSet(position[Index2D(ig,jg)].x(), position[Index2D(ig,jg)].y());
		}
		

        // GROUND COLLISION
        if ( position[Index2D(ig,jg)].y() <= (double) 0.05 ){
            hasItCollided = true;
            collisionNormal = Vector2d( 0, 1);
        }
        else if(position[Index2D(ig,jg)].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1, 0);
        }
        else if( position[Index2D(ig,jg)].y() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(0, -1.0);
        }
        else if( position [Index2D(ig,jg)].x() >= (double) 0.95 ){
            hasItCollided = true;
            collisionNormal = Vector2d(-1.0, 0);
        }



        if ( hasItCollided ){

            Vector2d velocityRel = velocity[Index2D(ig,jg)] - velocityCollision;

            double normalVelocity = velocityRel.dot(collisionNormal);

            if (normalVelocity < 0){

                Vector2d tangentVelocity = velocityRel - collisionNormal*normalVelocity;
                Vector2d velocityRelPrime;

                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity ){
                    velocityRelPrime = Vector2d (0, 0);
                }
                else{
                    velocityRelPrime.x() = tangentVelocity.x() + frictionCoeff*normalVelocity*tangentVelocity.x()/tangentVelocity.norm();
                    velocityRelPrime.y() = tangentVelocity.y() + frictionCoeff*normalVelocity*tangentVelocity.y()/tangentVelocity.norm();
            }

                velocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

            } // if normalVelocity < 0

        } // if hasItCollided




        //            cout << "Grid Node ( " << ig << " , " <<  jg << " , " << kg << " ) = " <<  hasItCollided << endl;

    } // for massList

}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
//====================================================================================================//
// GRID_HAT_OBSTACLE
//====================================================================================================//
/////////////////////////////////////////////////////////////////////////////////////////////////////////




double GRID_HAT_OBSTACLE::LevelSet(const double x, const double y){
	if (x >= 0.4 && x <= 0.6){
		return 0.5 - fabs(x-0.5)-y;
	}	
	else{
		return 1;
	}
}
Vector2d GRID_HAT_OBSTACLE::GradientLevelSet(const double x, const double y){
	Vector2d v;
	if ( x >= 0.4 && x <= 0.5 ){
		v = Vector2d( -1, 1 );
	}
	else if (x > 0.5 && x <= 0.6 ){
		v = Vector2d( 1, 1 );
	}
	return v/v.norm();
}

void GRID_HAT_OBSTACLE::GridCollision( double frictionCoeff ) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d ( 0, 0);

    for (int i = 0; i < massList.size(); i++){

        bool hasItCollided = false;

        int ig, jg;
        ig = massList[i].x();
        jg = massList[i].y();
		
		
		
		//////////////////////////////////////////
		// COLLISION WITH CYLINDER
		//////////////////////////////////////////
		double phi = LevelSet( position[Index2D(ig,jg)].x(), position[Index2D(ig,jg)].y() );
		if ( phi <= 0 ){
			// cout << "Collision with Cylinder GRID" << endl;
			hasItCollided = true;
			collisionNormal = GradientLevelSet(position[Index2D(ig,jg)].x(), position[Index2D(ig,jg)].y());
		}
		

        // GROUND COLLISION
        if ( position[Index2D(ig,jg)].y() <= (double) 0.05 ){
            hasItCollided = true;
            collisionNormal = Vector2d( 0, 1);
        }
        else if(position[Index2D(ig,jg)].x() <= (double) 0.05){
            hasItCollided = true;
            collisionNormal = Vector2d(1, 0);
        }
        else if( position[Index2D(ig,jg)].y() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(0, -1.0);
        }
        else if( position [Index2D(ig,jg)].x() >= (double) 0.95 ){
            hasItCollided = true;
            collisionNormal = Vector2d(-1.0, 0);
        }



        if ( hasItCollided ){

            Vector2d velocityRel = velocity[Index2D(ig,jg)] - velocityCollision;

            double normalVelocity = velocityRel.dot(collisionNormal);

            if (normalVelocity < 0){

                Vector2d tangentVelocity = velocityRel - collisionNormal*normalVelocity;
                Vector2d velocityRelPrime;

                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity ){
                    velocityRelPrime = Vector2d (0, 0);
                }
                else{
                    velocityRelPrime.x() = tangentVelocity.x() + frictionCoeff*normalVelocity*tangentVelocity.x()/tangentVelocity.norm();
                    velocityRelPrime.y() = tangentVelocity.y() + frictionCoeff*normalVelocity*tangentVelocity.y()/tangentVelocity.norm();
            }

                velocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

            } // if normalVelocity < 0

        } // if hasItCollided




        //            cout << "Grid Node ( " << ig << " , " <<  jg << " , " << kg << " ) = " <<  hasItCollided << endl;

    } // for massList

}