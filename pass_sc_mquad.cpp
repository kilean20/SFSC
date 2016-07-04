#include "pass_sc_mquad.h"
#include "pass_sc_drift.h"
#include "pass_sc_gaussian.h"
#include "global.h"
#include "statvec.h"

#include <cmath>
using namespace std;

void statVecMap_mQuadKick(STATvec & sigma, double KL)
{
    STATvec sigma2=sigma;
   // npx+npz=1
    for (unsigned nx=0;nx<3;nx++){
        for (unsigned nz=0;nz<3-nx;nz++){
            sigma2(nx,1,nz,0)-= KL*sigma(nx+1,0,nz,0);
            sigma2(nx,0,nz,1)+= KL*sigma(nx,0,nz+1,0);
        }
    }
    // npx+npz=2
    for (unsigned nx=0;nx<2;nx++){
        for (unsigned nz=0;nz<2-nx;nz++){
            sigma2(nx,2,nz,0)+= KL*( -2.0*sigma(nx+1,1,nz,0) +KL*sigma(nx+2,0,nz,0) );
            sigma2(nx,0,nz,2)+= KL*( 2.0*sigma(nx,0,nz+1,1) +KL*sigma(nx,0,nz+2,0) );
            sigma2(nx,1,nz,1)-= KL*( sigma(nx+1,0,nz,1) -sigma(nx,1,nz+1,0)
                                      +KL*sigma(nx+1,0,nz+1,0) );
        }
    }
    // npx+npz=3
    for (unsigned nx=0;nx<1;nx++){
        for (unsigned nz=0;nz<1-nx;nz++){
            sigma2(nx,3,nz,0)+= KL*( -3.0*sigma(nx+1,2,nz,0)
                                      +KL*(3.0*sigma(nx+2,1,nz,0)
                                            -KL*sigma(nx+3,0,nz,0)));
            sigma2(nx,0,nz,3)+= KL*( 3.0*sigma(nx,0,nz+1,2)
                                      +KL*(3.0*sigma(nx,0,nz+2,1)
                                            +KL*sigma(nx,0,nz+3,0)));
            sigma2(nx,2,nz,1)+= KL*( -2.0*sigma(nx+1,1,nz,1)+sigma(nx,2,nz+1,0)
                                        +KL*(sigma(nx+2,0,nz,1)-2.0*sigma(nx+1,1,nz+1,0)
                                              +KL*sigma(nx+2,0,nz+1,0)));
            sigma2(nx,1,nz,2)+= KL*( 2.0*sigma(nx,1,nz+1,1)-sigma(nx+1,0,nz,2)
                                        +KL*(sigma(nx,1,nz+2,0)-2.0*sigma(nx+1,0,nz+1,1)
                                              -KL*sigma(nx+1,0,nz+2,0)));
        }
    }
    // nx+nz=4
    sigma2(0,4,0,0)+= KL*(-4.0*sigma(1,3,0,0)
                         +KL*(6.0*sigma(2,2,0,0)
                             +KL*(-4.0*sigma(3,1,0,0)
                                 +KL*sigma(4,0,0,0) )));
    sigma2(0,0,0,4)+= KL*(4.0*sigma(0,0,1,3)
                           +KL*(6.0*sigma(0,0,2,2)
                               +KL*(4.0*sigma(0,0,3,1)
                                   +KL*sigma(0,0,4,0) )));
    sigma2(0,3,0,1)+= KL*(sigma(0,3,1,0)-3.0*sigma(1,2,0,1)
                         +KL*(-3.0*sigma(1,2,1,0)+3.0*sigma(2,1,0,1)
                             +KL*(3.0*sigma(2,1,1,0)-sigma(3,0,0,1)
                                 -KL*sigma(3,0,1,0) )));
    sigma2(0,1,0,3)+= KL*(-sigma(1,0,0,3)+3.0*sigma(0,1,1,2)
                         +KL*(-3.0*sigma(1,0,1,2)+3.0*sigma(0,1,2,1)
                             +KL*( -3.0*sigma(1,0,2,1) +sigma(0,1,3,0)
                                 -KL*sigma(1,0,3,0) )));
    sigma2(0,2,0,2)+= 2.0*KL*(sigma(0,2,1,1)-sigma(1,1,0,2)
                         +KL*(0.5*sigma(0,2,2,0)-2.0*sigma(1,1,1,1)+0.5*sigma(2,0,0,2)
                             +KL*(-sigma(1,1,2,0)+sigma(2,0,1,1)
                                 +0.5*KL*sigma(2,0,2,0) )));
    sigma=sigma2;
}



void sc_mQuadPass(STATvec & sigma, double  L, double K1, double  Ksc, unsigned Nint){

    const double ds = 0.5*L/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            statVecMap_mQuadKick(sigma, R4[0]*ds*K1);
        }
        //drift
        statVecMap_Drift(sigma,R4[0]*ds);
        //SC kick
        statVecMap_SC_gaussian(sigma, 2.0*R4[0]*ds*Ksc);
        //drift
        statVecMap_Drift(sigma, R4[0]*ds);

        for(unsigned r=1; r<3; r++)
        {
            //kick
            statVecMap_mQuadKick(sigma, (R4[r]+R4[r-1])*ds*K1);
            //drift
            statVecMap_Drift(sigma,R4[r]*ds);
            //SC kick
            statVecMap_SC_gaussian(sigma, 2.0*R4[r]*ds*Ksc);
            //drift
            statVecMap_Drift(sigma,R4[r]*ds);
        }
        //kick
        if(i==Nint-1){
            statVecMap_mQuadKick(sigma, R4[2]*ds*K1);
        }else{
            statVecMap_mQuadKick(sigma, 2.0*R4[2]*ds*K1);
        }
    }

//    //drift
//    statVecMap_Drift(sigma,0.5*ds);
//    for (unsigned i=0;i<Nint-1;i++){
//        //quadkick
//        statVecMap_mQuadKick(sigma, 0.5*ds*K1);
//        //SC kick
//        statVecMap_SC_gaussian(sigma, ds*Ksc);
//        //quadkick
//        statVecMap_mQuadKick(sigma, 0.5*ds*K1);
//        //drift
//        statVecMap_Drift(sigma,ds);
//    }
//    //quadkick
//    statVecMap_mQuadKick(sigma, 0.5*ds*K1);
//    //SC kick
//    statVecMap_SC_gaussian(sigma, ds*Ksc);
//    //quadkick
//    statVecMap_mQuadKick(sigma, 0.5*ds*K1);
//    //drift
//    statVecMap_Drift(sigma,0.5*ds);
}
