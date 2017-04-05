//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#include <vector>
#include "../Eigen/Dense"
#include <iostream>

//#include "GRID.h"

#include "../include/SIMULATION_PARAMETERS.h"


using namespace Eigen;
using namespace std;


////////////////////////////////////////////////////////////////////////
//**********************************************************************
// SETTING THE PARAMETERS CLASS... ALL OTHER PARAMETERS ARE INHERITED
//**********************************************************************
////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////
// SETTING THE TIME STEP / CFL CONDITION
//////////////////////////////////////////////////////////////

void SIMULATION_PARAMETERS::SetSimulationName(string simulationName) {
    m_name = simulationName;
}

void SIMULATION_PARAMETERS::SetDt(double timeStep) {
    m_dt = timeStep;
}

void SIMULATION_PARAMETERS::SetDt(const double h, const double timeToFrame, const vector<Vector2d> &particleVelocity) {
    VectorXd normParticleVelocity(particleVelocity.size());

    for (int p = 0; p < particleVelocity.size(); p++){
        normParticleVelocity[p] = particleVelocity[p].norm();
    }

    double maxVp = normParticleVelocity.maxCoeff();
    // cout << "Max Velocity = " << maxVp << endl;

    if (maxVp < 1e-5) {
        m_dt = 0.0001;
    }
    else {
        m_dt = m_CFL*h/maxVp;
		// if (m_dt > 1e-4){
		// 	m_dt = 1e-3;
		// }
		if (m_dt > timeToFrame){
			m_dt = timeToFrame;
		}
    }

}



//////////////////////////////////////////////////////////////////
// THE GET FUNCTIONS
//////////////////////////////////////////////////////////////////

string SIMULATION_PARAMETERS::GetSimulationName() {
    return m_name;
}

double SIMULATION_PARAMETERS::GetTimeToFrame(){
	return m_timeToFrame;
}

double SIMULATION_PARAMETERS::GetDt() {
    return m_dt;
}

double SIMULATION_PARAMETERS::GetFinalTime() {
    return m_finalTime;
}

double SIMULATION_PARAMETERS::GetAlpha()  {
    return m_alpha;
}

double SIMULATION_PARAMETERS::GetBeta() {
    return m_beta;
}

double SIMULATION_PARAMETERS::GetCriticalCompression() {
    return m_criticalCompression;
}

double SIMULATION_PARAMETERS::GetCriticalStretch() {
    return m_criticalStretch;
}

double SIMULATION_PARAMETERS::GetHardeningCoeff() {
    return m_hardeningCoeff;
}

double SIMULATION_PARAMETERS::GetInitialDensity() {
    return m_initialDensity;
}

double SIMULATION_PARAMETERS::GetYoungsMod() {
    return m_youngsMod;
}

double SIMULATION_PARAMETERS::GetPoisson() {
    return m_poissonRatio;
}

double SIMULATION_PARAMETERS::GetLambda() {
    return m_lambda0;
}

double SIMULATION_PARAMETERS::GetMu() {
    return m_mu0;
}

double SIMULATION_PARAMETERS::GetFrictionCoeff() {
    return m_frictionCoeff;
}

bool SIMULATION_PARAMETERS::SolveImplicit() {
    return m_implicitSolve;
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// DEFAULT PARAMETERS
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void DEFAULT_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.5;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.6;

    m_name = string("default");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}


void DEFAULT_IMPLICIT_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.5;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 1; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = true;


    m_CFL = 0.6;

    m_name = string("default_implicit");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// LOWER YOUNGS MODULUS PARAMETERS
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void LOWER_YOUNGS_MODULUS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.50;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 4.8e+4;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("lower_youngs_modulus");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}







void LOWER_YOUNGS_MODULUS_IMPLICIT::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.50;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 1; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 4.8e+4;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = true;


    m_CFL = 0.1;

    m_name = string("lower_youngs_modulus_implicit");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// LOWER HARDENING
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void LOWER_HARDENING::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.50;
	m_timeToFrame = 1.0/60.0;
    
	m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 5;
   m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("lower_hardening");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}





///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// LOWER_CRITICAL_COMPRESSION_PARAMETERS
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////








void LOWER_CRITICAL_COMPRESSION_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.5;
	m_timeToFrame = 1.0/60.0;
	
    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 1.9e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("lower_critical_compression");

    gravity = Vector2d (0,-9.8);
    usePlasticity = true;
}








///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// LOWER_CRITICAL_STRETCH_PARAMETERS
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////





void LOWER_CRITICAL_STRETCH_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.5;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


   	m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 5.0e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("lower_critical_stretch");

    gravity = Vector2d (0,-9.8);
    usePlasticity = true;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS
//
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.50;
	m_timeToFrame = 1.0/60.0;
	
    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 1.9e-2;
    m_criticalStretch = 5.0e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 1.4e+5;
	// m_youngsMod = 2;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("lower_critical_compression_stretch");


    gravity = Vector2d (0,-9.8);
    usePlasticity = true;


}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void HYPERELASTICITY::SetDefaultParameters() {
    m_dt = 0.0001;
    m_finalTime = 2.5;
	m_timeToFrame = 1.0/60.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.7;

    // Material Coefficients
    m_criticalCompression = 1.9e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 1.0e+2;
    m_youngsMod = 4.0e+4;
    m_poissonRatio = 0.2;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.1;

    m_name = string("hyperelastic");
    gravity = Vector2d (0,-9.8);

    usePlasticity = false;

}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



void SCALED_DEFAULT_PARAMETERS::SetDefaultParameters() {
    m_dt = 0.00001;
    m_finalTime = 2.00;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.5;

    // Material Coefficients
    m_criticalCompression = 2.5e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 3;
    m_initialDensity = 4.0e+2;
    // m_youngsMod = 1.4e+5;
	m_youngsMod = 2;
    m_poissonRatio = 0.3;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.3;

    m_name = string("Default_Scaled");


    gravity = Vector2d (0,-1);
    usePlasticity = true;


}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




void SCALED_HYPERELASTICITY::SetDefaultParameters() {
    m_dt = 0.00001;
    m_finalTime = 2.0;

    m_alpha = 0.95; // Grid to Particle Velocity Transfer PIC/FLIP

    m_beta = 0; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    m_frictionCoeff = 0.5;

    // Material Coefficients
    m_criticalCompression = 1.9e-2;
    m_criticalStretch = 7.5e-3;
    m_hardeningCoeff = 10;
    m_initialDensity = 4.0e+2;
    // m_youngsMod = 2.4e+5;
	m_youngsMod = 2;
    m_poissonRatio = 0.3;

    // Initial Lame Parameters
    m_lambda0 = m_youngsMod*m_poissonRatio/( ( 1+m_poissonRatio)*(1 - 2*m_poissonRatio) );
    m_mu0 = m_youngsMod/( 2*(1 + m_poissonRatio) );
    m_implicitSolve = false;


    m_CFL = 0.3;

    m_name = string("Hyperelastic_Scaled");
    gravity = Vector2d (0,-1);

    usePlasticity = false;

}