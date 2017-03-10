//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#ifndef MPM_LAB_2D_ELASTOPLASTIC_H
#define MPM_LAB_2D_ELASTOPLASTIC_H

#include <vector>
#include "../Eigen/Dense"
#include <iostream>

using namespace std;
using namespace Eigen;

// COMPUTING THE CONSTITUTIVE MODEL
// PSI is the strain energy


class ELASTOPLASTIC{

public:
    vector<double> ComputeMu(bool usePlasticity, double mu0, double hardeningCoeff, vector<double>& JP);

    vector<double> ComputeLambda(bool usePlasticity, double lambda0, double hardeningCoeff, vector<double>& JP);

    vector<double> Psi(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<double>& JElastic,
                       vector<double>& JPlastic, vector<Matrix2d>& elasticDeformationGradient, vector<Matrix2d>& R,
                       vector<Matrix2d>& S);

    vector<Matrix2d> dPsidF( bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<double>& JElastic,
                             vector<double>& JPlastic, vector<Matrix2d>& elasticDeformationGradient,
                             vector<Matrix2d>& R, vector<Matrix2d>& S);

    void CauchyStress(bool usePlasticity, double mu0, double lambda0, double hardeningCoeff, vector<Matrix2d>& sigma, vector<double>& JElastic,
                      vector<double>& JPlastic, vector<Matrix2d>& elasticDeformationGradient,
                      vector<Matrix2d>& R, vector<Matrix2d>& S);

    void ComputePolarDecomposition( vector<Matrix2d>& F, vector<Matrix2d>& R, vector<Matrix2d>& S );
};


#endif //MPM_LAB_2D_ELASTOPLASTIC_H
