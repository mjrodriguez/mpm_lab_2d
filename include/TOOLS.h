//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#ifndef MPM_LAB_2D_TOOLS_H
#define MPM_LAB_2D_TOOLS_H

#include <vector>
#include "../Eigen/Dense"
#include "PARTICLES.h"



using namespace std;
using namespace Eigen;


class TOOLS{



public:

    double Sum( vector<double>& input);
    int Index3D(const int N, const int i, const int j, const int k);
    int Index2D(const int N, const int i, const int j);
//    void ComputePolarDecomposition(PARTICLES& Particle);
    double MaxNormValue( vector<Vector3d>& V);
    double MinNormValue( vector<Vector3d>& V);
    double MaxNormValue( vector<Matrix3d>& V);
    double MinNormValue( vector<Matrix3d>& V);


    double MaxNormValue( vector<Vector2d>& V);
    double MinNormValue( vector<Vector2d>& V);
    double MaxNormValue( vector<Matrix2d>& V);
    double MinNormValue( vector<Matrix2d>& V);

};




#endif //MPM_LAB_2D_TOOLS_H
