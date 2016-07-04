//#include <iostream>
#include <fstream>
#include <cmath>
#include "pass_sc_mquad.h"
#include "pass_sc_equad.h"
#include "pass_sc_drift.h"
#include "pass_rk4_sc_drift.h"
#include "statvec.h"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <include/armadillo>
#else
    #include <armadillo>
#endif

using namespace std;
using namespace arma;

int main()
{
    ofstream Ksc0;//,Ksc1,Ksc2;
    Ksc0.open ("emittance_Ksc0.data");
    //Ksc1.open ("Ksc1e9_62resonance.data");
    //Ksc2.open ("Ksc1e8_62resonance.data");

    //cout.precision(15);

    //Define the ring
    const double ld=1.0;
    const double lq=0.5;
    const double kf=1.369, kd=1.427;

    double Ksc=1e-14;
    const double emittance = 0.115E-5,  betx=3.50146, alfx=-1.32962;
    const double gamx=(1+alfx*alfx)/betx;

    unsigned nTurn=20000;
    STATvec sigma;
    sigma(2,0,0,0)=emittance*betx;
    sigma(1,1,0,0)=emittance*alfx;
    sigma(0,2,0,0)=emittance*gamx;
    sigma(0,0,2,0)=emittance*betx;
    sigma(0,0,1,1)=-emittance*alfx;
    sigma(0,0,0,2)=emittance*gamx;

    sigma(4,0,0,0)=3.0*emittance*emittance*betx*betx;
    sigma(3,1,0,0)=3.0*emittance*emittance*betx*alfx;
    sigma(2,2,0,0)=emittance*emittance*(betx*gamx+2.0*alfx*alfx);
    sigma(1,3,0,0)=3.0*emittance*emittance*gamx*alfx;
    sigma(0,4,0,0)=3.0*emittance*emittance*gamx*gamx;

    sigma(0,0,4,0)=3.0*emittance*emittance*betx*betx;
    sigma(0,0,3,1)=-3.0*emittance*emittance*betx*alfx;
    sigma(0,0,2,2)=emittance*emittance*(betx*gamx+2.0*alfx*alfx);
    sigma(0,0,1,3)=-3.0*emittance*emittance*gamx*alfx;
    sigma(0,0,0,4)=3.0*emittance*emittance*gamx*gamx;


    for(unsigned i=0;i<nTurn;i++){
        for(unsigned j=0;j<18;j++){
            sc_Drift_Pass(sigma, 0.5*ld, Ksc, 50);
            sc_mQuadPass(sigma, lq, kf, Ksc, 100);
            sc_Drift_Pass(sigma, ld, Ksc, 100);
            sc_mQuadPass(sigma, lq, -kd, Ksc, 100);
            sc_Drift_Pass(sigma, 0.5*ld, Ksc, 50);
        }
        Ksc0 << sigma(2,0,0,0)*sigma(0,2,0,0)-sigma(1,1,0,0)*sigma(1,1,0,0) << endl;
        //Ksc0 << sigma(2,0,0,0) << "   " << sigma(0,2,0,0) << "   " << sigma(0,0,2,0) << "   " << sigma(0,0,0,2) << endl;
    }

//    //Ksc = 1e-9;
//    sigma.zeros();
//    sigma(2,0,0,0)=1e-4;
//    sigma(0,2,0,0)=1e-4;
//    sigma(0,0,2,0)=1e-5;
//    sigma(0,0,0,2)=1e-5;
//    Ksc=1e-9;
//    for(unsigned i=0;i<nTurn;i++){
//        for(unsigned j=0;j<50;j++){
//            sc_Drift_Pass(sigma, 0.5*lDrift_, Ksc, 6);
//            sc_mQuadPass(sigma, lQuad, quadK, Ksc, 20);
//            sc_Drift_Pass(sigma, lDrift_, Ksc, 12);
//            sc_mQuadPass(sigma, lQuad, -quadK, Ksc, 20);
//            sc_Drift_Pass(sigma, 0.5*lDrift_, Ksc, 6);
//        }
//        Ksc1 << sigma(2,0,0,0) << "   " << sigma(0,2,0,0) << "   " << sigma(0,0,2,0) << "   " << sigma(0,0,0,2) << endl;
//    }

//    //Ksc = 1e-8;
//    sigma.zeros();
//    sigma(2,0,0,0)=1e-4;
//    sigma(0,2,0,0)=1e-4;
//    sigma(0,0,2,0)=1e-5;
//    sigma(0,0,0,2)=1e-5;
//    Ksc = 1e-8;
//    for(unsigned i=0;i<nTurn;i++){
//        for(unsigned j=0;j<50;j++){
//            sc_Drift_Pass(sigma, 0.5*lDrift_, Ksc, 6);
//            sc_mQuadPass(sigma, lQuad, quadK, Ksc, 20);
//            sc_Drift_Pass(sigma, lDrift_, Ksc, 12);
//            sc_mQuadPass(sigma, lQuad, -quadK, Ksc, 20);
//            sc_Drift_Pass(sigma, 0.5*lDrift_, Ksc, 6);
//        }
//        Ksc2 << sigma(2,0,0,0) << "   " << sigma(0,2,0,0) << "   " << sigma(0,0,2,0) << "   " << sigma(0,0,0,2) << endl;
//    }

    Ksc0.close();//Ksc1.close();Ksc2.close();

    return 0;
}
