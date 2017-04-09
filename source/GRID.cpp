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

        //
        // force[i].x() += mass[ Index2D(ig, jg) ]*gravity.x();
        // force[i].y() += mass[ Index2D(ig, jg) ]*gravity.y();
		
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


void GRID::UpdateGravityVelocityGrid(const double timeStep, const Vector2d& gravity) {
    for (int i = 0; i < massList.size(); i++){
        int ig, jg;

        ig = massList[i].x();
        jg = massList[i].y();

		newVelocity[ Index2D(ig, jg) ] += timeStep*gravity;
		
    }
}





///////////////////////////////////////////////////////////////////





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
            A = newVelocity[ Index2D(ig, jg) ]*gradientWeight.transpose();


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









/////////////////////////////////////////////////////////////////////////////////////
// IMPLICIT TIME INTEGRATION MODULE
// SOLVES SYSTEM OF EQUATION DESCRIBED BY EQUATION (9) IN STOMAKHIN ET AL. 2013.
// ALSO DESCRIBED IN THESIS - SEMI IMPLICIT INTEGRATION SYSTEM
// 
// SUBMODULES:
// 		1. COMPUTE deltaF
//		2. COMPUTE delta_JF_iT
//		3. COMPUTE delta_R
//		4. COMPUTE Ap
//		
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
// PERFORM THIS AFTER COMPUTING THE EXPLICIT UPDATE
// EXPLICIT UPDATE SHOULD BE YOUR INITIAL GUESS
//////////////////////////////////////////////////////////
void GRID::ImplicitUpdateVelocityGrid(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic){
	
	// Setting the velocity grid nodes to temporary vector u, then we can compute the
	// hessian for an arbitrary vector u rather than specify for velocity and an 
	// arbitrary vector.
	vector<Vector2d> u;
	
	for (int i = 0; i < massList.size(); i++){
		int ig = massList[i].x();
		int jg = massList[i].y();
		
		u.push_back(newVelocity[ Index2D(ig,jg) ]);
	}
	
	
	// vector<Vector2d> V1;
	// Vector2d v(0,1);
	// V1.push_back(v);
	// V1.push_back(v);
	//
	// cout << V1[0].transpose() << endl;
	// cout << V1[1].transpose() << endl;
	// cout << "Inner Product^2 of V1 = " << InnerProduct(V1,V1) << endl;
	// cin.get();

	// void ConjugateResidual(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic);


	
	ConjugateResidual(usePlasticity, beta, dt, mu0, lambda0, hardeningCoeff, JE, JP, particleVolume, u, positionParticle, FE, R, S, Interpolation, Elastoplastic);
	
	for (int i = 0; i < massList.size(); i++){
		int ig = massList[i].x();
		int jg = massList[i].y();
		
		newVelocity[Index2D(ig,jg)] = u[i];
	}
	
	
}

void GRID::ConjugateResidual(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic){
	

	
	// Max Iterations
	int iter = 0;
	int maxIters = 30;
	double  alpha, y, test;
	// compute r -> residual
	vector<Vector2d> r, s, p, q;
	
	// letting Eu = r temporarily to save some space...
	r = ComputeHessianAction(usePlasticity, beta, dt, mu0, lambda0, hardeningCoeff, JE, JP, particleVolume, u, positionParticle, FE, R, S, Interpolation, Elastoplastic); 

	r = VectorSubtraction(u, r);
	
	for (int i = 0; i < u.size(); i++){
		r[i] = u[i] - r[i];
	}
	
	
	// for (int i = 0; i < u.size(); i++){
	// 	cout << u[i].transpose() << endl;
	// }
	// cin.get();
	// s = Ar
	s = ComputeHessianAction(usePlasticity, beta, dt, mu0, lambda0, hardeningCoeff, JE, JP, particleVolume, r, positionParticle, FE, R, S, Interpolation, Elastoplastic); 
	p = r;
	q = s;

	y = InnerProduct(r, s);
	alpha = y/ InnerProduct(q,q);
	test = sqrt( InnerProduct(r,r) );
	
	while ( test > sqrt(InnerProduct(u,u))*1e-5 && iter <= maxIters ){
			
		cout << "Iteration " << iter << " | " << "Norm = " << test << endl;
		double alpha_num = InnerProduct(r,s);
		double alpha_den = InnerProduct(q,q);
		alpha = alpha_num/alpha_den;
		double BETA_den = alpha_num; // Save for steps below;
		
		// UPDATE SOLUTION AND RESIDUAL
		for (int i = 0; i < u.size(); i++){
			u[i] += alpha*p[i];
			r[i] -= alpha*q[i];
		}
		
		s = ComputeHessianAction(usePlasticity, beta, dt, mu0, lambda0, hardeningCoeff, JE, JP, particleVolume, r, positionParticle, FE, R, S, Interpolation, Elastoplastic);
		
		double BETA_num = InnerProduct(r,s);
		double BETA = BETA_num/BETA_den;
		
		// UPDATE p AND q
		for (int i = 0; i < p.size(); i++){
			p[i] = r[i] + BETA*p[i];
			q[i] = s[i] + BETA*q[i];
		}
		
		test = sqrt( InnerProduct(r,r) );
		iter += 1;
		
		// vector<Vector2d> V = scalarVectorProduct(alpha,p);
		//
		// u = VectorAddition(u, V);
		// V.clear();
		//
		// V = scalarVectorProduct(alpha,q);
		// r = VectorSubtraction(r, V);
		// V.clear();
		// // Updating s = Ar with
		// // r_j+1 which was updated in the line above.
		// s = ComputeHessianAction(usePlasticity, beta, dt, mu0, lambda0, hardeningCoeff, JE, JP, particleVolume, r, positionParticle, FE, R, S, Interpolation, Elastoplastic);
		//
		// double B = InnerProduct(r, s)/y;
		// y = B*y;
		// V = scalarVectorProduct(B,p);
		// p = VectorAddition( r, V );
		// V.clear();
		//
		// V = scalarVectorProduct(B,q);
		// q = VectorAddition( s, V);
		// V.clear();
		//
		// test = fabs( alpha * InnerProduct(r,r) );
		// iter += 1;
	}

	
}

// NOT REALLY BEING USED....
vector<Vector2d> GRID::scalarVectorProduct(double a, vector<Vector2d>& p){
	vector<Vector2d> temp;
	for (int i = 0; i < p.size(); i++){
		temp.push_back( a*p[i] );
	}
	return temp;
}
vector<Vector2d> GRID::VectorSubtraction(vector<Vector2d>& p, vector<Vector2d>& q){
	vector<Vector2d> temp;
	for (int i = 0; i < p.size(); i++){
		temp.push_back( p[i]-q[i] );
	}
	return temp;
}
vector<Vector2d> GRID::VectorAddition(vector<Vector2d>& p, vector<Vector2d>& q){
	vector<Vector2d> temp;
	for (int i = 0; i < p.size(); i++){
		temp.push_back( p[i]+q[i] );
	}
	return temp;
}





double GRID::InnerProduct(vector<Vector2d>& a, vector<Vector2d>& b){
	double sum = 0;
	for (int i = 0; i < a.size(); i++){		
		sum += a[i].x()*b[i].x() + a[i].y()*b[i].y(); // A[i]*B[i] + A[2*i+1]*B[2*i+1];
	}
	return sum;
	
}


vector<Vector2d> GRID::ComputeHessianAction(const bool usePlasticity, const double beta, const double dt, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<double>& particleVolume, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& FE, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic){
	vector<Vector2d> Au, deltaForce;
	vector<Matrix2d> Ap;
	// ComputeAp(const double timeStep, const bool usePlasticity, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<Vector2d>& positionParticle, const vector<Vector2d>& velocityGrid, const vector<Vector2d>& u, const vector<Matrix2d>& F, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic)
	
	// 
	Ap = ComputeAp(dt, usePlasticity, mu0, lambda0, hardeningCoeff, JE, JP, positionParticle, u, FE, R, S, Interpolation, Elastoplastic);
	deltaForce = ComputeDeltaForce(particleVolume, positionParticle, FE, Ap, Interpolation );
	for (int i = 0; i < deltaForce.size(); i++){
		double ig = massList[i].x();
		double jg = massList[i].y();
		double one_over_mass = 1.0/mass[ Index2D(ig,jg) ];
		
		Au.push_back( u[i] - beta*dt*deltaForce[i]*one_over_mass );
	}
	
	return Au;
	
}






vector<Vector2d> GRID::ComputeDeltaForce(const vector<double>& particleVolume, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& F, const vector<Matrix2d>& Ap, INTERPOLATION& Interpolation  ){
    vector<Vector2d> delta_force;
	Vector2d DELTA_FORCE;
	for (int i = 0; i < massList.size(); i++){
        
        
		DELTA_FORCE = Vector2d::Zero();
		
        int ig, jg;
        ig = massList[i].x();
        jg = massList[i].y();
			
		for (int p = 0; p < F.size(); p++){
            Vector2d gradientWeight;
			gradientWeight = Vector2d::Zero();
			
			gradientWeight = Interpolation.GradientWeight(m_h, positionParticle[p], ig, jg);
			DELTA_FORCE += particleVolume[p]*Ap[p]*F[p].transpose()*gradientWeight; 
			
		}
		
		delta_force.push_back(DELTA_FORCE);
		
	}
	
	return delta_force;

}




////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// vector<Matrix3d> GRID::ComputeDeltaF(const double timeStep, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& elasticDeformationGradient, INTERPOLATION& Interpolation)
// Input: timeStep ===> time step size, dt.
//        positionParticle ==> particle position vector<Vector2d>
//		  elasticDeformationGradient ==> deformation gradient FE
//		  INTERPOLATION ==> Interpolation class
// output: deltaF ==> vector<Matrix2d>
//
// Details: This ComputeDeltaF is used to compute hessian action Ev where v is the velocity
//
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////




// vector<Matrix2d> GRID::ComputeDeltaF(const double timeStep, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& elasticDeformationGradient, INTERPOLATION& Interpolation){
// 	vector<Matrix2d> deltaF;
//
//
// 	for (int p = 0; p < elasticDeformationGradient.size(); p++){
//         Vector2d gradientWeight;
//         Matrix2d deltaFp;
//
// 		for (int i = 0; i < massList.size(); i++){
//             int ig, jg;
//             ig = massList[i].x();
//             jg = massList[i].y();
//
// 			gradientWeight = Interpolation.GradientWeight(m_h, positionParticle[p], ig, jg);
// 			deltaFp += timeStep*velocity[ Index2D(ig, jg) ]*gradientWeight.transpose()*elasticDeformationGradient[p];
//
// 		}
//
// 		deltaF.push_back(deltaFp);
//
// 	}
//
// 	return deltaF;
//
// }



////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// vector<Matrix3d> GRID::ComputeDeltaF(const double timeStep, const vector<Vector2d>& u, const vector<Vector2d>& positionParticle, const vector<Matrix2d>& elasticDeformationGradient, INTERPOLATION& Interpolation)
// Input: timeStep ===> time step size, dt.
//        u ==> arbitrary vector u of size = massList.size()
//        positionParticle ==> particle position vector<Vector2d>
//		  elasticDeformationGradient ==> deformation gradient FE
//		  INTERPOLATION ==> Interpolation class
// output: deltaF ==> vector<Matrix2d>
//
// Details: This ComputeDeltaF is used to compute arbitrary hessian action Eu.
//
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////



vector<Matrix2d> GRID::ComputeDeltaF(const double timeStep, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, vector<double>& JE_hat, const vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& F_hat, vector<Matrix2d>& R_hat, vector<Matrix2d>& S_hat , INTERPOLATION& Interpolation){
	vector<Matrix2d> deltaF;
	
	// For computing Fhat
	Matrix2d velocity_Gradient, A;
	velocity_Gradient = Matrix2d::Zero();
	A = Matrix2d::Zero();
	
	for (int p = 0; p < elasticDeformationGradient.size(); p++){
        Matrix2d deltaFp;
		
		for (int i = 0; i < u.size(); i++){
            int ig, jg;
            ig = massList[i].x();
            jg = massList[i].y();
			
			Vector2d gradientWeight = Interpolation.GradientWeight(m_h, positionParticle[p], ig, jg);
			deltaFp += timeStep*u[i]*gradientWeight.transpose()*elasticDeformationGradient[p];
			
        	A = newVelocity[ Index2D(ig, jg) ]*gradientWeight.transpose();
		
        	velocity_Gradient += A  ;
			
		}
		
		deltaF.push_back(deltaFp);
// const double dt, const Matrix2d& velocityGradient, vector<double>& JE_hat, const vector<Vector2d>& particlePosition, const Matrix2d& F, vector<Matrix2d>& F_hat, vector<Matrix2d>&  R_hat, vector<Matrix2d>&  S_hat, INTERPOLATION& Interpolation
		ComputeFHat(timeStep, velocity_Gradient, JE_hat, positionParticle, elasticDeformationGradient[p], F_hat, R_hat, S_hat, Interpolation); 
		
	}
	
	return deltaF;
	
}






vector<Matrix2d> GRID::ComputeAp(const double timeStep, const bool usePlasticity, const double mu0, const double lambda0, const double hardeningCoeff, vector<double>& JE, vector<double>& JP, const vector<Vector2d>& positionParticle, vector<Vector2d>& u, const vector<Matrix2d>& F, const vector<Matrix2d>& R, const vector<Matrix2d>& S, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic){
    vector<Matrix2d> Ap, deltaF, deltaR, deltaJFiT, F_hat, R_hat, S_hat;
    vector<double> mu, lambda, JE_hat;
	
	// const double timeStep, vector<Vector2d>& u, const vector<Vector2d>& positionParticle, vector<double>& JE_hat, const vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& F_hat, vector<Matrix2d>& R_hat, vector<Matrix2d>& S_hat , INTERPOLATION& Interpolation
    deltaF = ComputeDeltaF(timeStep, u, positionParticle, JE_hat, F, F_hat, R_hat, S_hat, Interpolation);
	Elastoplastic.ComputePolarDecomposition(F_hat, R_hat, S_hat);
	
    deltaR = ComputeDeltaR(deltaF, S_hat, R_hat);
    deltaJFiT = ComputeDeltaJFiT(F_hat, deltaF);
    mu = Elastoplastic.ComputeMu(usePlasticity, mu0, hardeningCoeff, JP);
    lambda = Elastoplastic.ComputeLambda(usePlasticity, lambda0, hardeningCoeff, JP);


	
    for (int p = 0; p < F.size(); p ++){
        Matrix2d B = JE_hat[p]*F[p].transpose().inverse();
        double JFiT_dd_deltaF = DoubleDotProduct( B, deltaF[p]);
        Matrix2d TEMP = 2*mu[p]*deltaF[p] - 2*mu[p]*deltaR[p] + lambda[p]*JE_hat[p]*F_hat[p].transpose().inverse()*JFiT_dd_deltaF + lambda[p]*(JE_hat[p]-1)*deltaJFiT[p];
        Ap.push_back(TEMP);
    }




    return Ap;


}



void GRID::ComputeFHat(const double dt, const Matrix2d& velocityGradient, vector<double>& JE_hat, const vector<Vector2d>& particlePosition, const Matrix2d& F, vector<Matrix2d>& F_hat, vector<Matrix2d>&  R_hat, vector<Matrix2d>&  S_hat, INTERPOLATION& Interpolation){
	Matrix2d F_hat_temp = F + dt*velocityGradient*F;
	double JE_hat_temp = F_hat_temp.determinant();
	F_hat.push_back(F_hat_temp);
	JE_hat.push_back( JE_hat_temp );
	R_hat.push_back(Matrix2d::Zero());
	S_hat.push_back(Matrix2d::Zero());	
}




//
// // OLD VERSION
// void GRID::ComputeFHat(double dt, vector<double>& JE_hat, const vector<Vector2d>& particlePosition, const vector<Matrix2d>& F, vector<Matrix2d>& F_hat, vector<Matrix2d>&  R_hat, vector<Matrix2d>&  S_hat, INTERPOLATION& Interpolation, ELASTOPLASTIC& Elastoplastic){
//
// 	for (int p = 0; p < particlePosition.size(); p++){
// 		Vector2d gradientWeight;
// 		Matrix2d velocity_Gradient, A;
// 		velocity_Gradient = Matrix2d::Zero();
// 		A = Matrix2d::Zero();
//
// 		for ( int i = 0; i < massList.size(); i++ ) {
// 			int ig, jg;
// 			ig = massList[i].x();
// 			jg = massList[i].y();
// 			gradientWeight = Interpolation.GradientWeight(m_h, particlePosition[p], ig, jg);
//         	A = newVelocity[ Index2D(ig, jg) ]*gradientWeight.transpose();
//
//         	velocity_Gradient += A  ;
//     	} // END WEIGHT INTERPOLATION TRANSFER
// 		F_hat.push_back(F[p] + dt*velocity_Gradient*F[p]);
// 		JE_hat.push_back( F_hat[p].determinant() );
// 		R_hat.push_back(Matrix2d::Zero());
// 		S_hat.push_back(Matrix2d::Zero());
//
// 	}
//
// 	Elastoplastic.ComputePolarDecomposition(F_hat, R_hat, S_hat);
//
// }
//



////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// vector<Matrix3d> GRID::ComputeDeltaR(vector<Matrix3d> &deltaF, vector<Matrix3d> &S, vector<Matrix3d> &R)
// Input: deltaF ===> differential of deformation gradient, 3x3 matrix for all particles
//        S ==> polar decomposition, 3x3 matrix for all particles
//        R ==> polar decomposition matrix, 3x3 matrix for all particles
// output: deltaR ==> vector<Matrix3d>
//
// Details: Check thesis Implicit time integration section
//
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////

vector<Matrix2d> GRID::ComputeDeltaR(const vector<Matrix2d> &deltaF, const vector<Matrix2d> &S, const vector<Matrix2d> &R) {

    vector<Matrix2d> deltaR;
    for (int p = 0; p < S.size(); p++){

        Matrix2d B = R[p].transpose()*deltaF[p] - deltaF[p].transpose()*R[p];
        // b << B(0,1), B(0,2), B(1,2);
        // A << S[p](0,0)+S[p](1,1), S[p](2,1), -S[p](0,2),    S[p](1,2), S[p](0,0)+S[p](2,2), S[p](0,1),     -S[p](0,2), S[p](1,0), S[p](1,1)+S[p](2,2);
        // x =  A.inverse()*b; // A.colPivHouseholderQr().solve(b);
		
		double a = B(0,1)/(S[p](0,0) + S[p](1,1));
		
		Matrix2d RTdR;
       	RTdR << 0 , a, -a, 0;

        deltaR.push_back( R[p]*RTdR );
        RTdR = Matrix2d::Zero();

    }


    return deltaR;

}


















////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// vector<Matrix2d> GRID::ComputeDeltaJFiT(vector<Matrix2d> &deltaF, vector<Matrix2d> &S, vector<Matrix2d> &R)
// Input: deltaF ===> differential of deformation gradient, 3x3 matrix for all particles
//        S ==> polar decomposition, 2x2 matrix for all particles
//        R ==> polar decomposition matrix, 2x2 matrix for all particles
// output: deltaR ==> vector<Matrix3d>
//
// Details: Check thesis Implicit time integration section
//
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////



vector<Matrix2d> GRID::ComputeDeltaJFiT(const vector<Matrix2d>& F, const vector<Matrix2d>& deltaF){

        vector<Matrix2d> deltaJFiT;
        vector<Matrix2d> A; // 4x1 vector of 2x2 matrices
        for (int p = 0; p < F.size(); p++){
            A = Compute_dJFiT( F[p] );
            deltaJFiT.push_back(DoubleDotProduct(A,deltaF[p]));
            for (int i = 0; i < 4; i++) {
                A[i] = Matrix2d::Zero();
            }
        }

    return deltaJFiT;

}




////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// Compute_dJFiT(Matrix2d& Fe) -> derivative of J*inv(F)^T with respect components of F
// Input: Fe ===> Elastic deformation gradient
// output: Fourth order tensor (Vector of 4x1 2x2 matrices / Total of 16 components)
// 
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////

vector<Matrix2d> GRID::Compute_dJFiT(const Matrix2d& F){
	vector<Matrix2d> dJFiT;
	
	Matrix2d temp;
	
	// dJFiT(0)
	temp << 0, 0, 0, 1;
	dJFiT.push_back( temp );
	
	// dJFiT(1)
	temp << 0, 0, -1, 0;
	dJFiT.push_back( temp );
	
	temp << 0, -1, 0, 0;
	dJFiT.push_back( temp );
	
	temp << 1, 0, 0, 0;
	dJFiT.push_back( temp );
	
	
	return dJFiT;
								
	
}





////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// DoubleDotProduct(vector<Matrix3d>& C, Matrix3d& B)
// Input: 	C ===> Fourth Order Tensor in 4x1 of 2x2 matrices format.
//			B ===> second order tensor 2x2 matrix.
// output: 	A ===> second order tensor 2x2 matrix.
// 
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////

Matrix2d GRID::DoubleDotProduct(const vector<Matrix2d>& C, const Matrix2d& B){
	Matrix2d matrixProduct;
	Matrix2d A;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			matrixProduct = C[2*i + j]*B;
			A(i,j) = matrixProduct.sum();
		}
	}

    return A;
}







////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************************
// DoubleDotProduct(Matrix2d& C, Matrix2d& B)
// Input: 	C ===> second Order Tensor 2x2 matrix.
//			B ===> second order tensor 2x2 matrix.
// output: 	A ===> a double number in reals.
//
// Author: Martin Rodriguez Cruz
// Date: February 16, 2017
// ************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////

double GRID::DoubleDotProduct(const Matrix2d& C, const Matrix2d& B){
    Matrix2d matrixProduct;
    double A;

    matrixProduct = C*B;
    A = matrixProduct.sum();


    return A;
}
















////////////////////////////////////////////////////////////////////////////////////////////
//=======================================================================================//
// GRID_NO_OBSTACLE
//=======================================================================================//
///////////////////////////////////////////////////////////////////////////////////////////






void GRID_NO_OBSTACLE::GridCollision( double frictionCoeff ) {
    Vector2d collisionNormal;
    Vector2d velocityCollision = Vector2d ( 0, 0);

    for (int i = 0; i < massList.size(); i++){

        bool hasItCollided = false;
		bool stickyCollision = true;
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
			stickyCollision = true;
        }
        else if( position[Index2D(ig,jg)].y() >= (double) 0.95){
            hasItCollided = true;
            collisionNormal = Vector2d(0, -1.0);
			stickyCollision = true;
        }
        else if( position [Index2D(ig,jg)].x() >= (double) 0.95 ){
            hasItCollided = true;
            collisionNormal = Vector2d(-1.0, 0);
			stickyCollision = true;
        }



        if ( hasItCollided ){

            Vector2d velocityRel = newVelocity[Index2D(ig,jg)] - velocityCollision;

            double normalVelocity = velocityRel.dot(collisionNormal);

            if (normalVelocity < 0){

                Vector2d tangentVelocity = velocityRel - collisionNormal*normalVelocity;
                Vector2d velocityRelPrime;

                if (tangentVelocity.norm() <= -frictionCoeff*normalVelocity || stickyCollision ){
                    velocityRelPrime = Vector2d (0, 0);
                }
                else{
                    velocityRelPrime.x() = tangentVelocity.x() + frictionCoeff*normalVelocity*tangentVelocity.x()/tangentVelocity.norm();
                    velocityRelPrime.y() = tangentVelocity.y() + frictionCoeff*normalVelocity*tangentVelocity.y()/tangentVelocity.norm();
            }

                newVelocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

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
	return (x-0.5)*(x-0.5) + (y-0.3)*(y-0.3) - 0.1*0.1;
}
Vector2d GRID_CYLINDER_OBSTACLE::GradientLevelSet(const double x, const double y){
	Vector2d v = Vector2d( 2*(x-0.5), 2*(y-0.3) );
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

            Vector2d velocityRel = newVelocity[Index2D(ig,jg)] - velocityCollision;

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

                newVelocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

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
	if (x >= 0.4 && x <= 0.6 && y >= 0.5 && y <= 0.6){
		return -( 0.6 - fabs(x-0.5) - y );
	}	
	else{
		return (double) 1;
	}
}
Vector2d GRID_HAT_OBSTACLE::GradientLevelSet(const double x, const double y){
	Vector2d v;
	if ( x >= 0.4 && x <= 0.5 && y >= 0.5 && y <= 0.6 ){
		v = Vector2d( -1, 1 );
	}
	else if (x > 0.5 && x <= 0.6 && y >= 0.5 && y <= 0.6 ){
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

            Vector2d velocityRel = newVelocity[Index2D(ig,jg)] - velocityCollision;

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

                newVelocity[ Index2D(ig,jg) ] = velocityRelPrime + velocityCollision;

            } // if normalVelocity < 0

        } // if hasItCollided




        //            cout << "Grid Node ( " << ig << " , " <<  jg << " , " << kg << " ) = " <<  hasItCollided << endl;

    } // for massList

}
