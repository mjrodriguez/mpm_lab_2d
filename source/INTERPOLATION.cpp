//
// Created by Martin Rodriguez Cruz on 9/19/16.
//


#include "../Eigen/Dense"
#include "../include/INTERPOLATION.h"

using namespace Eigen;


double INTERPOLATION::N(const double x) {

    if (fabs(x) >= 0.0 && fabs(x) < 1.0) {
        return 0.5 * fabs(x) * fabs(x) * fabs(x) - x * x + 2.0 / 3.0;
    } else if (fabs(x) >= 1.0 && fabs(x) < 2.0) {
        return -1.0 / 6.0 * fabs(x) * fabs(x) * fabs(x) + x * x - 2 * fabs(x) + 4.0 / 3.0;
    } else {
        return 0.0;
    }

}


double INTERPOLATION::GradientN(const double x) {
    if ( fabs(x) >= 0.0 && fabs(x)  < 1.0 ){
        return 1.5*fabs(x)*x - 2.0*x;
    }
    else if ( fabs(x) >= 1.0 && fabs(x) < 2.0 ){
        return -0.5*x*fabs(x) + 2.0*x - 2.0*x/fabs(x);
    }
    else {
        return 0.0;
    }
	
	/*if (x >= -2.0 && x < -1.0){
		return 0.5* x * x + 2.0*x + 2.0;
	}
	else if (x >= -1.0 && x < 0.0 ) {
		return -1.5*x*x - 2.0*x;
	}
	else if (x >= 0.0 && x < 1.0){
		return 1.5*x*x - 2.0*x;
	}
	else if (x >= 1.0 && x <= 2.0){
		return -0.5*x*x + 2.0*x - 2.0;
	}
	else {
		return 0.0; 
	}*/
	
	
}


double INTERPOLATION::Weight(const double h, Vector3d &particlePosition, const double i, const double j, const double k) {
    return N( ( particlePosition.x() - i*h )/h ) * N( (particlePosition.y()-j*h)/h ) * N( (particlePosition.z()-k*h)/h );
}



Vector3d INTERPOLATION::GradientWeight(const double h, Vector3d &particlePosition, const double i, const double j, const double k) {
    Vector3d DW;
    DW[0] = GradientN( (particlePosition.x()-i*h)/h ) * N( (particlePosition.y()-j*h)/h ) * N( (particlePosition.z()-k*h)/h );
    DW[1] = N( (particlePosition.x()-i*h)/h ) * GradientN( (particlePosition.y()-j*h)/h ) * N( (particlePosition.z()-k*h)/h );
    DW[2] = N( (particlePosition.x()-i*h)/h ) * N( (particlePosition.y()-j*h)/h ) * GradientN( (particlePosition.z()-k*h)/h );

    return DW;

}







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MPM 2D - INTERPOLATION
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////








double INTERPOLATION::Weight(const double h, const Vector2d& particlePosition, const double i, const double j){
    return N( ( particlePosition.x() - i*h )/h ) * N( (particlePosition.y()-j*h)/h );
}

Vector2d INTERPOLATION::GradientWeight( const double h, const Vector2d &particlePosition, const double i, const double j){
    Vector2d DW;
    DW[0] = GradientN( (particlePosition.x()-i*h)/h ) * N( (particlePosition.y()-j*h)/h ) ;
    DW[1] = N( (particlePosition.x()-i*h)/h ) * GradientN( (particlePosition.y()-j*h)/h ) ;

    return DW;
}

