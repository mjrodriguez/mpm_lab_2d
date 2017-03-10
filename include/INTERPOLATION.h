//
// Created by Martin Rodriguez Cruz on 9/18/16.
//

#ifndef MPM_V4_0_INTERPOLATION_H
#define MPM_V4_0_INTERPOLATION_H

#include "../Eigen/Dense"

using namespace Eigen;

class INTERPOLATION{

public:

    double N(const double x);
    double GradientN(const double x);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // MPM 3D - INTERPOLATION
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double Weight( const double h, Vector3d& particlePosition, const double i, const double j, const double k);
    Vector3d GradientWeight( const double h, Vector3d& particlePosition, const double i, const double j, const double k);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // MPM 2D - INTERPOLATION
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double Weight( const double h, const Vector2d& particlePosition, const double i, const double j );
    Vector2d GradientWeight( const double h, const Vector2d &particlePosition, const double i, const double j);


};

#endif //MPM_V4_0_INTERPOLATION_H
