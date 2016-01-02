#ifndef GLOBAL_H
#define GLOBAL_H
#include <math.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <include/armadillo>
#else
    #include <armadillo>
#endif
using namespace arma;
//=========================================================================
//                           Global  constants 
//=========================================================================
//const double PI      = 3.14159265;
//const double light_speed =  299792458;
//--------------------------spin-parameters--------------------------------
const double MDM = 0.00115965218;  // anormalous magnetic momentum
const double EDM = 7.3E-17;  // electric momentum g-factor

//-------------------------------------------------------------------------
const double C1 = 1.0 / (2.0 - pow(2.0, 1.0/3.0));
const double C0 = 1.0 - 2.0*C1;

//-------------------------state-vector-index------------------------------
enum x_index {x_=0, px_=1, s_= 2, ps_=3, z_=4, pz_=5, vt_=6, dE_=7, Sx_=8, Ss_=9};

//-------------------------state-vector-index------------------------------
enum type_index {DRIFT_=0, eBEND_=1, eQUAD_= 2};

//-----------relativistic parameters for electron magic energy-------------
const double GAMMA=29.38243572993826;
const double BETA=0.9994206777202302;
const double BETA2=BETA*BETA;
const double BETAGAMMA=GAMMA*BETA; // beta_0 gamma_0
const double BETAGAMMA2=BETAGAMMA*BETAGAMMA;
const double magicENERGY = 15.014392631143503; //MeV
//-----------
#endif
                