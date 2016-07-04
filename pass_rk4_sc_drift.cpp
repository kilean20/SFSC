#include "pass_rk4_sc_drift.h"
#include "global.h"
#include "statvec.h"

#include <cmath>
using namespace std;

inline double
dx_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom){
    double temp=0.0;
    nx-=1;
    temp += sigma(nx,npx+1,nz,npz);

    if (nx+npx+nz+npz < 3)
        temp += //0.5*( sigma(nx,npx+3,nz,npz) +sigma(nx,npx+1,nz,npz+2) )
                + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx+1,nz,npz) + sigma(nx,npx+1,nz+2,npz)/r );
    return temp;
}

inline double
dz_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom){
    double temp=0.0;
    nz-=1;
    temp += sigma(nx,npx,nz,npz+1);
    if (nx+npx+nz+npz < 3)
        temp += //0.5*( sigma(nx,npx+2,nz,npz+1) +sigma(nx,npx,nz,npz+3) )
                + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx,nz,npz+1) + sigma(nx,npx,nz+2,npz+1)/r );
    return temp;
}

inline double
dpx_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom){
    double temp=0.0;
    npx-=1;
    temp += Ksc*denom*sigma(nx+1,npx,nz,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom*GAMMA2*(sigma(nx+1,npx+2,nz,npz)+sigma(nx+1,npx,nz,npz+2))
                - Ksc*denom*denom*( (2.0+r)/6*sigma(nx+3,npx,nz,npz) + 0.5/r*sigma(nx+1,npx,nz+2,npz) ) ;

    return temp;
}

inline double
dpz_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom){
    double temp=0.0;
    npz-=1;
    temp += Ksc*denom/r*sigma(nx,npx,nz+1,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom/r*GAMMA2*(sigma(nx,npx+2,nz+1,npz)+sigma(nx,npx,nz+1,npz+2))
                - Ksc*denom*denom/r*( 0.5*sigma(nx+2,npx,nz+1,npz) + (1.0+2.0*r)/6/r/r*sigma(nx,npx,nz+3,npz) ) ;

    return temp;
}

void sc_Drift_Pass_RK4(STATvec & sigma, double  L, double  Ksc, unsigned Nint){

    const double ds = L/(double)Nint;

    double r, denom;
    STATvec sigma1,sigma2,sigma3,sigma4,sigmaTemp;

    for (unsigned i=0;i<Nint;i++){
        denom =1.0/(sqrt(sigma(2,0,0,0))*(sqrt(sigma(2,0,0,0))+sqrt(sigma(0,0,2,0))));
        r=sqrt(sigma(0,0,2,0)/sigma(2,0,0,0));

        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma1(nx,npx,nz,npz)=0.5*ds*(nx*dx_Drift_RK4(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + npx*dpx_Drift_RK4(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + nz*dz_Drift_RK4(nx,npx,nz,npz,sigma,r,Ksc,denom)
                                                          + npz*dpz_Drift_RK4(nx,npx,nz,npz,sigma,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp= sigma+sigma1;
        denom =1.0/(sqrt(sigmaTemp(2,0,0,0))*(sqrt(sigmaTemp(2,0,0,0))+sqrt(sigmaTemp(0,0,2,0))));
        r=sqrt(sigmaTemp(0,0,2,0)/sigmaTemp(2,0,0,0));

        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma2(nx,npx,nz,npz)=0.5*ds*(nx*dx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npx*dpx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + nz*dz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npz*dpz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma2;
        denom =1.0/(sqrt(sigmaTemp(2,0,0,0))*(sqrt(sigmaTemp(2,0,0,0))+sqrt(sigmaTemp(0,0,2,0))));
        r=sqrt(sigmaTemp(0,0,2,0)/sigmaTemp(2,0,0,0));
        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma3(nx,npx,nz,npz)=ds*(nx*dx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npx*dpx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + nz*dz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                          + npz*dpz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma3;
        denom =1.0/(sqrt(sigmaTemp(2,0,0,0))*(sqrt(sigmaTemp(2,0,0,0))+sqrt(sigmaTemp(0,0,2,0))));
        r=sqrt(sigmaTemp(0,0,2,0)/sigmaTemp(2,0,0,0));
        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma4(nx,npx,nz,npz)=ds*(nx*dx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + npx*dpx_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + nz*dz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom)
                                                      + npz*dpz_Drift_RK4(nx,npx,nz,npz,sigmaTemp,r,Ksc,denom));
                        }
                    }
                }
            }
        }
        for (unsigned i=0;i<69;i++)    sigma(i)+=(sigma1(i)+2.0*sigma2(i)+sigma3(i)+0.5*sigma4(i))/3.0;
    }
}

