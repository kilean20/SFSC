#include <iostream>
#include "statvec.h"
#include "pass_sc_mquad.h"
#include "pass_sc_equad.h"
#include "pass_sc_drift.h"

using namespace std;


int main()
{
    const double lq=1.0;
    const double kf=1.369, kd=1.427;
    cout.precision(7);

    double Ksc=0.0;
    const double emittance = 0.115E-4,  betx=2.33135, alfx=-1.29714;
    const double gamx=(1+alfx*alfx)/betx;

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


    sc_mQuadPass(sigma, lq, kf, Ksc, 50);
    //sc_eQuadPass(sigma, 0.75, 0.0, 10);
    //sc_Drift_Pass(sigma, 0.75, 0, 10);
    cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;

    sigma.zeros();
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
    sigma(2,0,0,0)=emittance*betx;
    sigma(1,1,0,0)=emittance*alfx;
    sigma(0,2,0,0)=emittance*gamx;
    sigma(0,0,2,0)=emittance*betx;
    sigma(0,0,1,1)=-emittance*alfx;
    sigma(0,0,0,2)=emittance*gamx;
    sc_eQuadPass(sigma, lq, kf, Ksc, 50);
    cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;

    return 0;
}
