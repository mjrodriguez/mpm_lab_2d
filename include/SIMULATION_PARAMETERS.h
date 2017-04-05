//
// Created by Martin Rodriguez Cruz on 3/7/17.
//

#ifndef MPM_LAB_2D_SIMULATION_PARAMETERS_H
#define MPM_LAB_2D_SIMULATION_PARAMETERS_H

#include <vector>
#include "../Eigen/Dense"
#include <iostream>

using namespace Eigen;
using namespace std;



// ALL SIMULATIONS WILL USE THIS LIST OF PARAMETERS. TO CREATE A NEW LIST OF SIMULATION PARAMETERS THEN INHERIT "class SIMULATION_PARAMETERS"

class SIMULATION_PARAMETERS{

protected:
    double m_dt;
    double m_finalTime;
	double m_timeToFrame;
    double m_alpha;  // Grid to Particle Velocity Transfer PIC/FLIP

    double m_beta; // Explicit Integration = 0
    // Trapezoidal Integration = 1/2
    // Implicit Integration = 1


    double m_frictionCoeff;

    // Material Coefficients
    double m_criticalCompression;
    double m_criticalStretch;
    double m_hardeningCoeff;
    double m_initialDensity;
    double m_youngsMod;
    double m_poissonRatio;

    // Initial Lame Parameters
    double m_lambda0;
    double m_mu0;
    bool m_implicitSolve;

    double m_CFL;
    string m_name;


public:

    Vector2d gravity;
    bool usePlasticity;

    virtual void SetDefaultParameters() = 0;
    void SetDt(double timeStep);
    void SetDt(const double h, const double timeToFrame, const vector<Vector2d>& particleVelocity);
    void SetSimulationName(string simulationName);

    /////////////////////////////////////////////////////
    // FUNCTIONS TO GET PARAMETERS
    ////////////////////////////////////////////////////
    double GetDt();
	double GetTimeToFrame();
    double GetFinalTime();
    double GetAlpha();
    double GetBeta();
    double GetCriticalCompression();
    double GetCriticalStretch();
    double GetHardeningCoeff();
    double GetInitialDensity();
    double GetYoungsMod();
    double GetPoisson();
    double GetLambda();
    double GetMu();
    double GetFrictionCoeff();
    bool SolveImplicit();
    string GetSimulationName();


};


//////////////////////////////////////////////////////////////////////
// DEFAULT PARAMETERS [FROM MPM PAPER]
//////////////////////////////////////////////////////////////////////

class DEFAULT_PARAMETERS: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};

class DEFAULT_IMPLICIT_PARAMETERS: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};

//////////////////////////////////////////////////////////////////////
// LOWER_YOUNGS_MODULUS [FROM MPM PAPER]
//////////////////////////////////////////////////////////////////////

class LOWER_YOUNGS_MODULUS: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};

class LOWER_YOUNGS_MODULUS_IMPLICIT: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};

//////////////////////////////////////////////////////////////////////
// LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS [FROM MPM PAPER]
//////////////////////////////////////////////////////////////////////

class LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};

//////////////////////////////////////////////////////////////////////
// LOWER_HARDENING [FROM MPM PAPER]
//////////////////////////////////////////////////////////////////////

class LOWER_HARDENING: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};




///////////////////////////////////////////////////////////////////////
// LOWER CRITICAL COMPRESSION PARAMETERS
///////////////////////////////////////////////////////////////////////


class LOWER_CRITICAL_COMPRESSION_PARAMETERS: public SIMULATION_PARAMETERS{
public:
    void SetDefaultParameters();
};


///////////////////////////////////////////////////////////////////////
// LOWER CRITICAL STRETCH PARAMETERS
///////////////////////////////////////////////////////////////////////


class LOWER_CRITICAL_STRETCH_PARAMETERS: public SIMULATION_PARAMETERS{
public:
    void SetDefaultParameters();
};





///////////////////////////////////////////////////////////////////////
// HYPERELASTICITY PARAMETERS
///////////////////////////////////////////////////////////////////////


class HYPERELASTICITY: public SIMULATION_PARAMETERS{
public:
    void SetDefaultParameters();
};



//////////////////////////////////////////////////////////////////////
// DEFAULT PARAMETERS [FROM MPM PAPER]
//////////////////////////////////////////////////////////////////////

class SCALED_DEFAULT_PARAMETERS: public SIMULATION_PARAMETERS {
public:
    void SetDefaultParameters();

};




///////////////////////////////////////////////////////////////////////
// SCALED HYPERELASTICITY PARAMETERS
///////////////////////////////////////////////////////////////////////


class SCALED_HYPERELASTICITY: public SIMULATION_PARAMETERS{
public:
    void SetDefaultParameters();
};



#endif //MPM_LAB_2D_SIMULATION_PARAMETERS_H
