//
// Created by Martin Rodriguez Cruz on 3/7/17.
//



#include <vector>
#include "../Eigen/Dense"
#include "../include/PARTICLES.h"
#include "../include/TOOLS.h"


using namespace std;
using namespace Eigen;

//typedef Matrix<double,Dynamic,1> VectorXd;

double TOOLS::Sum(vector<double> &input) {
    double temp = 0;

    for (int i = 0; i < input.size(); i++){
        temp += input[i];
    }

    return temp;
}



int TOOLS::Index3D(const int N, const int i, const int j, const int k) {
    return i*N*N + j*N + k;
}


int TOOLS::Index2D(const int N, const int i, const int j){
    return i*N + j;
}


double TOOLS::MaxNormValue(vector<Vector3d> &V) {
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.maxCoeff();
}

double TOOLS::MaxNormValue(vector<Matrix3d> &V) {
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.maxCoeff();
}

double TOOLS::MinNormValue(vector<Matrix3d>& V){
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.minCoeff();
}

double TOOLS::MinNormValue(vector<Vector3d>& V){
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.minCoeff();
}






/////////////////////////////////////////////////////////



double TOOLS::MaxNormValue(vector<Vector2d> &V) {
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.maxCoeff();
}

double TOOLS::MaxNormValue(vector<Matrix2d> &V) {
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.maxCoeff();
}

double TOOLS::MinNormValue(vector<Matrix2d>& V){
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.minCoeff();
}

double TOOLS::MinNormValue(vector<Vector2d>& V){
    VectorXd normV(V.size());
    for (int i = 0; i < V.size(); i++){
        normV[i] = V[i].norm();
    }

    return normV.minCoeff();
}