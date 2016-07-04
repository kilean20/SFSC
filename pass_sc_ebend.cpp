#include "pass_sc_ebend.h"
#include "global.h"
#include "statvec.h"

#include <cmath>
using namespace std;

double
dx_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom){
    double temp=0.0;
    nx-=1;
    temp += sigma(nx,npx+1,nz,npz);

    if (nx+npx+nz+npz < 4)
        temp += 0.5*iRho*sigma(nx+1,npx+1,nz,npz);
    if (nx+npx+nz+npz < 3)
        temp += 0.5*( sigma(nx,npx+3,nz,npz) +sigma(nx,npx+1,nz,npz+2) )
              + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx+1,nz,npz) + sigma(nx,npx+1,nz+2,npz)/r ) // SC part
              + iRho*iRho*(2.0-0.5*BETA2)*sigma(nx+2,npx+1,nz,npz);
    if (nx+npx+nz+npz < 2)
        temp += 2.0*iRho*(sigma(nx+1,npx+3,nz,npz) +sigma(nx+1,npx+1,nz,npz+2))
                + iRho*iRho*iRho*(7.0/3.0-1.5*BETA2)*sigma(nx+3,npx+1,nz,npz);

    return temp;
}

double
dz_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom){
    double temp=0.0;
    nz-=1;
    temp += sigma(nx,npx,nz,npz+1);

    if (nx+npx+nz+npz < 4)
        temp += 0.5*iRho*sigma(nx+1,npx,nz,npz+1);
    if (nx+npx+nz+npz < 3)
        temp += 0.5*( sigma(nx,npx+2,nz,npz+1) +sigma(nx,npx,nz,npz+3) )
                + 0.5*Ksc*denom*GAMMA2* ( sigma(nx+2,npx,nz,npz+1) + sigma(nx,npx,nz+2,npz+1)/r )
                + iRho*iRho*(2.0-0.5*BETA2)*sigma(nx+2,npx,nz,npz+1);
    if (nx+npx+nz+npz < 2)
        temp += 2.0*iRho*(sigma(nx+1,npx+2,nz,npz+1) +sigma(nx+1,npx,nz,npz+3))
                + iRho*iRho*iRho*(7.0/3.0-1.5*BETA2)*sigma(nx+3,npx,nz,npz+1);

    return temp;
}

double
dpx_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom){
    double temp=0.0;
    npx-=1;
    temp += (Ksc*denom - iRho*iRho*(2.0-BETA2))*sigma(nx+1,npx,nz,npz);

    if (nx+npx+nz+npz < 4)
        temp += -iRho*(sigma(nx,npx+2,nz,npz) +sigma(nx,npx,nz,npz+2))
                + iRho*iRho*iRho*(1.5*BETA2-1.0)*sigma(nx+2,npx,nz,npz);
    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom*GAMMA2*(sigma(nx+1,npx+2,nz,npz)+sigma(nx+1,npx,nz,npz+2))
                - Ksc*denom*denom*( (2.0+r)/6.0*sigma(nx+3,npx,nz,npz) + 0.5/r*sigma(nx+1,npx,nz+2,npz) )
                - iRho*iRho*( (2.0-0.5*BETA2)*( sigma(nx+1,npx+2,nz,npz) +sigma(nx+1,npx,nz,npz+2) )
                              + iRho*iRho/6*(10.0-11.0*BETA2+3.0*BETA2*BETA2)*sigma(nx+3,npx,nz,npz) );
    if (nx+npx+nz+npz < 2)
        temp -= iRho*(  iRho*iRho*iRho*iRho*(22.0-40.0*BETA2+15.0*BETA2*BETA2)/12.0*sigma(nx+4,npx,nz,npz)
                       +0.25*iRho*iRho*(14.0-9.0*BETA2)*(sigma(nx+2,npx+2,nz,npz) +sigma(nx+2,npx,nz,npz+2))
                       +0.5*sigma(nx,npx+4,nz,npz) +sigma(nx,npx+2,nz,npz+2) +0.5*sigma(nx,npx,nz,npz+4));

    return temp;
}

double
dpz_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom){
    double temp=0.0;
    npz-=1;
    temp += (Ksc*denom/r)*sigma(nx,npx,nz+1,npz);

    if (nx+npx+nz+npz < 3)
        temp += 0.5*Ksc*denom/r*GAMMA2*(sigma(nx,npx+2,nz+1,npz)+sigma(nx,npx,nz+1,npz+2))
                - Ksc*denom*denom/r*( 0.5*sigma(nx+2,npx,nz+1,npz) + (1.0+2.0*r)/6/r/r*sigma(nx,npx,nz+3,npz) );

    return temp;
}

void sc_eBendPass(STATvec & sigma, double  L, double Angle, double Ksc, unsigned Nint){

    const double iRho = L/Angle;
    const double ds = L/(double)Nint;

    double r, denom;
    STATvec sigma1,sigma2,sigma3,sigma4,sigmaTemp;

    for (unsigned i=0;i<Nint;i++){
        denom =1.0/(sigma(2,0,0,0)*(sigma(2,0,0,0)+sigma(0,0,2,0)));
        r=sigma(0,0,2,0)/sigma(2,0,0,0);

        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma1(nx,npx,nz,npz)=0.5*ds*(nx*dx_eBend(nx,npx,nz,npz,sigma,r,iRho,Ksc,denom)
                                                          + npx*dpx_eBend(nx,npx,nz,npz,sigma,r,iRho,Ksc,denom)
                                                          + nz*dz_eBend(nx,npx,nz,npz,sigma,r,iRho,Ksc,denom)
                                                          + npz*dpz_eBend(nx,npx,nz,npz,sigma,r,iRho,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma1;
        r=(sigma(0,0,2,0)+sigma1(0,0,2,0))/(sigma(2,0,0,0)+sigma1(2,0,0,0));
        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma2(nx,npx,nz,npz)=0.5*ds*(nx*dx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                          + npx*dpx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                          + nz*dz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                          + npz*dpz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma2;
        r=(sigma(0,0,2,0)+sigma2(0,0,2,0))/(sigma(2,0,0,0)+sigma2(2,0,0,0));
        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma3(nx,npx,nz,npz)=ds*(nx*dx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + npx*dpx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + nz*dz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + npz*dpz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom));
                        }
                    }
                }
            }
        }
        sigmaTemp=sigma+sigma3;
        r=(sigma(0,0,2,0)+sigma3(0,0,2,0))/(sigma(2,0,0,0)+sigma3(2,0,0,0));
        for (unsigned nx=0;nx<4;nx++){
            for (unsigned npx=0;npx<4;npx++){
                for (unsigned nz=0;nz<4;nz++){
                    for (unsigned npz=0;npz<4;npz++){
                        if(nx+npx+nz+npz<5){
                            sigma4(nx,npx,nz,npz)=ds*(nx*dx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + npx*dpx_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + nz*dz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom)
                                                      + npz*dpz_eBend(nx,npx,nz,npz,sigmaTemp,r,iRho,Ksc,denom));
                        }
                    }
                }
            }
        }
        for (unsigned i=0;i<70;i++)    sigma(i)+=(sigma1(i)+2.0*sigma2(i)+sigma3(i)+0.5*sigma4(i))/3.0;
    }
}

