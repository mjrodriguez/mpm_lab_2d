//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <vector>
#include "../Eigen/Dense"
#include <iostream>


#include "../include/ELASTOPLASTIC.h"

using namespace std;
using namespace Eigen;


vector<double> ELASTOPLASTIC::ComputeMu(bool usePlasticity, double mu0, double hardeningCoeff, vector<double> &JP) {
    vector<double> MU;
    if (usePlasticity  ){
        for (int p = 0; p < JP.size(); p++ ){
			double noise = rand()/RAND_MAX + 1;
            MU.push_back( (1+noise)*mu0*exp( hardeningCoeff*(1-JP[p]) ) );
        }
    }
    else{
        for (int p = 0; p < JP.size(); p++ ) {
            MU.push_back(mu0);
        }
    }
    return MU;
}


vector<double> ELASTOPLASTIC::ComputeLambda( bool usePlasticity, double lambda0, double hardeningCoeff, vector<double> &JP) {
    vector<double> LAMBDA;

    if ( usePlasticity ){
        for (int p = 0; p < JP.size(); p++ ){
			double noise = rand()/RAND_MAX + 1;
            LAMBDA.push_back( (1+noise)*lambda0*exp( hardeningCoeff*( 1-JP[p] ) ) );
        }
    }
    else {
        for (int p = 0; p < JP.size(); p++ ) {
            LAMBDA.push_back(lambda0);
        }
    }

    return LAMBDA;
}


void ELASTOPLASTIC::ComputePolarDecomposition(vector<Matrix2d>& F, vector<Matrix2d>& R, vector<Matrix2d>& S ) {
    for (int p = 0; p < F.size(); p++){
        JacobiSVD<Matrix2d> svd(F[p], ComputeFullU | ComputeFullV);
        R[p] = svd.matrixU()*svd.matrixV().transpose();
        S[p] = svd.matrixV()*svd.singularValues().asDiagonal()*svd.matrixV().transpose();
    }
}


void ELASTOPLASTIC::CauchyStress(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<Matrix2d>& sigma, vector<double>& JElastic, vector<double>& JPlastic, vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& R, vector<Matrix2d>& S) {

    /*vector<Matrix3d> diffPsi;

    for (int p = 0; p < JElastic.size(); p++){
        diffPsi.push_back( Matrix3d::Zero() );
    }


    diffPsi = dPsidF(usePlasticity, mu0, lambda0, hardeningCoeff, JElastic, JPlastic, elasticDeformationGradient, R, S);
     */

    Matrix2d Eye = Matrix2d::Identity();
	
	

    vector<double> mu = ComputeMu(usePlasticity, mu0, hardeningCoeff, JPlastic);
    vector<double> lambda = ComputeLambda(usePlasticity, lambda0, hardeningCoeff, JPlastic);
    ComputePolarDecomposition( elasticDeformationGradient, R, S );

    for (int p = 0; p < JElastic.size(); p++){
        // sigma[p] = diffPsi[p]*elasticDeformationGradient[p].transpose()/(JElastic[p]*JPlastic[p]);
        double one_over_J = 1/( JElastic[p] * JPlastic[p]);
        sigma[p] =  2* mu[p] * one_over_J* (elasticDeformationGradient[p] - R[p]) *elasticDeformationGradient[p].transpose() + lambda[p]*one_over_J*(JElastic[p] - 1)*JElastic[p]*Eye    ;
    }

    mu.erase(mu.begin(), mu.end());
    lambda.erase(lambda.begin(), lambda.end());
}