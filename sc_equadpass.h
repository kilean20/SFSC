#ifndef SC_EQUADPASS_H
#define SC_EQUADPASS_H

#include <cmath>
using namespace std;

#include "global.h"
#include "statvec.h"

double
dx_eQuad(size_t nx, size_t npx, size_t nz, size_t npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom){
    double temp=0.0;
    nx-=1;
    temp += sigma(nx,npx+1,nz,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*( sigma(nx,npx+3,nz,npz) +sigma(nx,npx+1,nz,npz+2) )
              + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx+1,nz,npz) + sigma(nx,npx+1,nz+2,npz)/r ) // SC part
              + 0.5*K1*( sigma(nx+2,npx+1,nz,npz) -sigma(nx,npx+1,nz+2,npz) );
    return temp;
}

double
dz_eQuad(size_t nx, size_t npx, size_t nz, size_t npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom){
    double temp=0.0;
    nz-=1;
    temp += sigma(nx,npx,nz,npz+1);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*( sigma(nx,npx+2,nz,npz+1) +sigma(nx,npx,nz,npz+3) )
                + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx,nz,npz+1) + sigma(nx,npx,nz+2,npz+1)/r )
                + 0.5*K1*( sigma(nx+2,npx,nz,npz+1) -sigma(nx,npx,nz+2,npz+1) );
    return temp;
}

double
dpx_eQuad(size_t nx, size_t npx, size_t nz, size_t npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom){
    double temp=0.0;
    npx-=1;
    temp += (Ksc*denom -K1)*sigma(nx+1,npx,nz,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom*GAMMA2*(sigma(nx+1,npx+2,nz,npz)+sigma(nx+1,npx,nz,npz+2))
                - Ksc*denom*denom*( (2.0+r)/6*sigma(nx+3,npx,nz,npz) + 0.5/r*sigma(nx+1,npx,nz+2,npz) )
                - 0.5*K1*( sigma(nx+1,npx+2,nz,npz) +sigma(nx+1,npx,nz,npz+2)
                           + K1/GAMMA2*( sigma(nx+3,npx,nz,npz) -sigma(nx+1,npx,nz+2,npz) ) );

    return temp;
}

double
dpz_eQuad(size_t nx, size_t npx, size_t nz, size_t npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom){
    double temp=0.0;
    npz-=1;
    temp += (Ksc*denom/r +K1)*sigma(nx,npx,nz+1,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom/r*GAMMA2*(sigma(nx,npx+2,nz+1,npz)+sigma(nx,npx,nz+1,npz+2))
                - Ksc*denom*denom/r*( 0.5*sigma(nx+2,npx,nz+1,npz) + (1.0+2.0*r)/6/r/r*sigma(nx,npx,nz+3,npz) )
                + 0.5*K1*( sigma(nx,npx+2,nz+1,npz) +sigma(nx,npx,nz+1,npz+2)
                           + K1/GAMMA2*( sigma(nx+2,npx,nz+1,npz) -sigma(nx,npx,nz+3,npz) ) );

    return temp;
}

void sc_eQuadPass(STATvec & sigma, double  L, double K1, double  Ksc, size_t Nint){

    const double ds = L/(double)Nint;

    double r, denom;
    STATvec sigma1,sigma2,sigma3,sigma4,sigmaTemp;

    for (size_t i=0;i<Nint;i++){
        denom =1.0/(sigma(2,0,0,0)*(sigma(2,0,0,0)+sigma(0,0,2,0)));
        r=sigma(0,0,2,0)/sigma(2,0,0,0);

        for (size_t nx=0;nx<4;nx++){
            for (size_t npx=0;npx<4;npx++){
                for (size_t nz=0;nz<4;nz++){
                    for (size_t npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma1(nx,npx,nz,npz)=0.5*ds*(nx*dx_eQuad(nx,npx,nz,npz,sigma,r,K1,Ksc,denom)
                                                          + npx*dpx_eQuad(nx,npx,nz,npz,sigma,r,K1,Ksc,denom)
                                                          + nz*dz_eQuad(nx,npx,nz,npz,sigma,r,K1,Ksc,denom)
                                                          + npz*dpz_eQuad(nx,npx,nz,npz,sigma,r,K1,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma1;
        r=(sigma(0,0,2,0)+sigma1(0,0,2,0))/(sigma(2,0,0,0)+sigma1(2,0,0,0));
        for (size_t nx=0;nx<4;nx++){
            for (size_t npx=0;npx<4;npx++){
                for (size_t nz=0;nz<4;nz++){
                    for (size_t npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma2(nx,npx,nz,npz)=0.5*ds*(nx*dx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                          + npx*dpx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                          + nz*dz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                          + npz*dpz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma2;
        r=(sigma(0,0,2,0)+sigma2(0,0,2,0))/(sigma(2,0,0,0)+sigma2(2,0,0,0));
        for (size_t nx=0;nx<4;nx++){
            for (size_t npx=0;npx<4;npx++){
                for (size_t nz=0;nz<4;nz++){
                    for (size_t npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma3(nx,npx,nz,npz)=ds*(nx*dx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + npx*dpx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + nz*dz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + npz*dpz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma3;
        r=(sigma(0,0,2,0)+sigma3(0,0,2,0))/(sigma(2,0,0,0)+sigma3(2,0,0,0));
        for (size_t nx=0;nx<4;nx++){
            for (size_t npx=0;npx<4;npx++){
                for (size_t nz=0;nz<4;nz++){
                    for (size_t npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma4(nx,npx,nz,npz)=ds*(nx*dx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + npx*dpx_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + nz*dz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom)
                                                      + npz*dpz_eQuad(nx,npx,nz,npz,sigmaTemp,r,K1,Ksc,denom));
                        }
                    }
                }
            }
        }
        for (size_t i=0;i<70;i++)    sigma(i)+=(sigma1(i)+2.0*sigma2(i)+sigma3(i)+0.5*sigma4(i))/3.0;
    }
}

#endif

