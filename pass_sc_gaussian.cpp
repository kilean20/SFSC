#include "pass_sc_gaussian.h"
#include "global.h"
#include "statvec.h"

//#include <iostream>
#include <cmath>
using namespace std;

void statVecMap_SC_gaussian(STATvec & sigma, double KL)
{
    STATvec sigma2=sigma;
    const double denom=1.0/(sqrt(sigma(2,0,0,0))*(sqrt(sigma(2,0,0,0))+sqrt(sigma(0,0,2,0))));
    const double r=sqrt(sigma(0,0,2,0)/sigma(2,0,0,0));
    const double KLn=KL*denom;
    // npx+npz=1
    for (unsigned nx=0;nx<3;nx++){
        for (unsigned nz=0;nz<3-nx;nz++){
            sigma2(nx,1,nz,0)+= KLn*sigma(nx+1,0,nz,0);
            sigma2(nx,0,nz,1)+= KLn/r*sigma(nx,0,nz+1,0);
            if (nx+nz < 2){
                sigma2(nx,1,nz,0)-= KLn*denom*( (2.0+r)/6.0*sigma(nx+3,0,nz,0) +0.5/r*sigma(nx+1,0,nz+2,0) );
                sigma2(nx,0,nz,1)-= KLn*denom/r*( (r+0.5)/3.0*sigma(nx,0,nz+3,0)/r/r
                                                     +0.5*sigma(nx+2,0,nz+1,0));
            }
        }
    }
    // npx+npz=2
    for (unsigned nx=0;nx<2;nx++){
        for (unsigned nz=0;nz<2-nx;nz++){
            sigma2(nx,2,nz,0)+= KLn*( 2.0*sigma(nx+1,1,nz,0) +KLn*sigma(nx+2,0,nz,0) );
            sigma2(nx,0,nz,2)+= KLn/r*( 2.0*sigma(nx,0,nz+1,1) +KLn*sigma(nx,0,nz+2,0)/r );
            sigma2(nx,1,nz,1)+= KLn/r*( r*sigma(nx+1,0,nz,1)+sigma(nx,1,nz+1,0)
                                      +KLn*sigma(nx+1,0,nz+1,0) );
        }
    }
    sigma2(0,2,0,0)-= KLn*denom/3.0*( (2.0+r)*sigma(3,1,0,0) +3.0*sigma(1,1,2,0)/r
                                +KLn*( (2.0+r)*sigma(4,0,0,0) +3.0*sigma(2,0,2,0)/r ));
    sigma2(0,0,0,2)-= KLn*denom/r/r/3.0*( (2.0*r+1.0)*sigma(0,0,3,1)/r +3.0*sigma(2,0,1,1)*r
                                +KLn*( (2.0+1.0/r)*sigma(0,0,4,0)/r +3.0*sigma(2,0,2,0) ));
    sigma2(0,1,0,1)-= KLn*denom/6.0*( (2.0+r)*sigma(3,0,0,1) +(1.0/r+2.0)*sigma(0,1,3,0)/r/r
                                      +3.0*(sigma(2,1,1,0)+sigma(1,0,2,1))/r
                                +KLn/r*( (5.0+r)*sigma(3,0,1,0) + (5.0+1.0/r)*sigma(1,0,3,0)/r ));
    // npx+npz=3
    for (unsigned nx=0;nx<1;nx++){
        for (unsigned nz=0;nz<1-nx;nz++){
            sigma2(nx,3,nz,0)+= KLn*( 3.0*sigma(nx+1,2,nz,0)
                                      +KLn*(3.0*sigma(nx+2,1,nz,0)
                                            +KLn*sigma(nx+3,0,nz,0)));
            sigma2(nx,0,nz,3)+= KLn/r*( 3.0*sigma(nx,0,nz+1,2)
                                      +KLn/r*(3.0*r*sigma(nx,0,nz+2,1)
                                            +KLn/r*sigma(nx,0,nz+3,0)));
            sigma2(nx,2,nz,1)+= KLn/r*( 2.0*r*sigma(nx+1,1,nz,1)+sigma(nx,2,nz+1,0)
                                        +KLn*(r*sigma(nx+2,0,nz,1)+2.0*sigma(nx+1,1,nz+1,0)
                                              +KLn*sigma(nx+2,0,nz+1,0)));
            sigma2(nx,1,nz,2)+= KLn/r*( 2.0*sigma(nx,1,nz+1,1)+r*sigma(nx+1,0,nz,2)
                                        +KLn/r*(sigma(nx,1,nz+2,0)+2.0*r*sigma(nx+1,0,nz+1,1)
                                              +KLn*sigma(nx+1,0,nz+2,0)));
        }
    }
    // nx+nz=4
    sigma2(0,4,0,0)+= KLn*(4.0*sigma(1,3,0,0)
                         +KLn*(6.0*sigma(2,2,0,0)
                             +KLn*(4.0*sigma(3,1,0,0)
                                 +KLn*sigma(4,0,0,0) )));
    sigma2(0,0,0,4)+= KLn/r*(4.0*sigma(0,0,1,3)
                           +KLn/r*(6.0*sigma(0,0,2,2)
                               +KLn/r*(4.0*sigma(0,0,3,1)
                                   +KLn/r*sigma(0,0,4,0) )));
    sigma2(0,3,0,1)+= KLn/r*(sigma(0,3,1,0)+3.0*r*sigma(1,2,0,1)
                         +KLn*(3.0*sigma(1,2,1,0)+3.0*r*sigma(2,1,0,1)
                             +KLn*(3.0*sigma(2,1,1,0)+r*sigma(3,0,0,1)
                                 +KLn*sigma(3,0,1,0) )));
    sigma2(0,1,0,3)+= KLn/r*(r*sigma(1,0,0,3)+3.0*sigma(0,1,1,2)
                         +KLn*(3.0*sigma(1,0,1,2)+3.0/r*sigma(0,1,2,1)
                             +KLn*(3.0/r*sigma(1,0,2,1)+sigma(0,1,3,0)/r/r
                                 +KLn*sigma(1,0,3,0)/r/r )));
    sigma2(0,2,0,2)+= 2.0*KLn/r*(sigma(0,2,1,1)+r*sigma(1,1,0,2)
                         +KLn*(0.5*sigma(0,2,2,0)/r+2.0*sigma(1,1,1,1)+0.5*r*sigma(2,0,0,2)
                             +KLn*(sigma(1,1,2,0)/r+sigma(2,0,1,1)
                                 +0.5*KLn*sigma(2,0,2,0)/r )));
    sigma=sigma2;
}

