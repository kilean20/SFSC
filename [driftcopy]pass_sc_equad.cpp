#include "pass_sc_equad.h"
#include "pass_sc_drift.h"
#include "global.h"
#include "statvec.h"

#include <cmath>
#include <iostream>
using namespace std;

void sc_eQuadPass(STATvec & sigma, double  L, double K1, double Ksc, short Nint){
    cout << "this is equad" << endl;
    Ksc += K1;
    const double ds = L/(double)Nint;

    double r, denom;
    STATvec sigma1,sigma2,sigma3,sigma4,sigmaTemp;

    for (short i=0;i<Nint;i++){
        denom =1.0/(sigma(2,0,0,0)*(sigma(2,0,0,0)+sigma(0,0,2,0)));
        r=sigma(0,0,2,0)/sigma(2,0,0,0);

        for (short nx=0;nx<4;nx++){
            for (short npx=0;npx<4;npx++){
                for (short nz=0;nz<4;nz++){
                    for (short npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma1(nx,npx,nz,npz)=0.5*ds*(nx*dx_Drift(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + npx*dpx_Drift(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + nz*dz_Drift(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + npz*dpz_Drift(nx,npx,nz,npz,sigma,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp= sigma+sigma1;
        r=(sigma(0,0,2,0)+sigma1(0,0,2,0))/(sigma(2,0,0,0)+sigma1(2,0,0,0));
        for (short nx=0;nx<4;nx++){
            for (short npx=0;npx<4;npx++){
                for (short nz=0;nz<4;nz++){
                    for (short npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma2(nx,npx,nz,npz)=0.5*ds*(nx*dx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npx*dpx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + nz*dz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npz*dpz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma2;
        r=(sigma(0,0,2,0)+sigma2(0,0,2,0))/(sigma(2,0,0,0)+sigma2(2,0,0,0));
        for (short nx=0;nx<4;nx++){
            for (short npx=0;npx<4;npx++){
                for (short nz=0;nz<4;nz++){
                    for (short npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma3(nx,npx,nz,npz)=ds*(nx*dx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npx*dpx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + nz*dz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npz*dpz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma3;
        r=(sigma(0,0,2,0)+sigma3(0,0,2,0))/(sigma(2,0,0,0)+sigma3(2,0,0,0));
        for (short nx=0;nx<4;nx++){
            for (short npx=0;npx<4;npx++){
                for (short nz=0;nz<4;nz++){
                    for (short npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma4(nx,npx,nz,npz)=ds*(nx*dx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + npx*dpx_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + nz*dz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + npz*dpz_Drift(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        for (short i=0;i<70;i++)    sigma(i)+=(sigma1(i)+2.0*sigma2(i)+sigma3(i)+0.5*sigma4(i))/3.0;
    }
}

