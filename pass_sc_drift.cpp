#include "pass_sc_drift.h"
#include "pass_sc_gaussian.h"
#include "global.h"
#include "statvec.h"

#include <cmath>
using namespace std;

inline void statVecMap_Drift(STATvec & sigma, double L)
{
    STATvec sigma2=sigma;
    // nx+nz=1
    for (unsigned npx=0;npx<3;npx++){
        for (unsigned npz=0;npz<3-npx;npz++){
            sigma2(1,npx,0,npz)+= L*sigma(0,npx+1,0,npz);
            sigma2(0,npx,1,npz)+= L*sigma(0,npx,0,npz+1);
//            if (npx+npz < 2)
//            {
//                sigma2(1,npx,0,npz)+= 0.5*L*( sigma(0,npx+3,0,npz) +sigma(0,npx+1,0,npz+2) );
//                sigma2(0,npx,1,npz)+= 0.5*L*( sigma(0,npx,0,npz+3) +sigma(0,npx+2,0,npz+1) );
//            }
        }
    }
    // nx+nz=2
    for (unsigned npx=0;npx<2;npx++){
        for (unsigned npz=0;npz<2-npx;npz++){
            sigma2(2,npx,0,npz)+= L*( 2.0*sigma(1,npx+1,0,npz)+L*sigma(0,npx+2,0,npz) );
            sigma2(0,npx,2,npz)+= L*( 2.0*sigma(0,npx,1,npz+1)+L*sigma(0,npx,0,npz+2) );
            sigma2(1,npx,1,npz)+= L*( sigma(1,npx,0,npz+1)+sigma(0,npx+1,1,npz)+L*sigma(0,npx+1,0,npz+1) );

        }
    }
//    sigma2(2,0,0,0)+= L*(sigma(1,3,0,0)+sigma(1,1,0,2)+L*( sigma(0,4,0,0)+sigma(0,2,0,2) ));
//    sigma2(0,0,2,0)+= L*(sigma(0,0,1,3)+sigma(0,2,1,1)+L*( sigma(0,0,0,4)+sigma(0,2,0,2) ));
//    sigma2(1,0,1,0)+= 0.5*L*(sigma(1,2,0,1)+sigma(1,0,0,3)+sigma(0,3,1,0)+sigma(0,1,1,2)
//                             +2.0*L*( sigma(0,3,0,1)+sigma(0,1,0,3) ));
    // nx+nz=3
    for (unsigned npx=0;npx<1;npx++){
        for (unsigned npz=0;npz<1-npx;npz++){
            sigma2(3,npx,0,npz)+= L*( 3.0*sigma(2,npx+1,0,npz) +L*(3.0*sigma(1,npx+2,0,npz)
                                                                   +L*sigma(0,npx+3,0,npz) ));
            sigma2(0,npx,3,npz)+= L*( 3.0*sigma(0,npx,2,npz+1) +L*(3.0*sigma(0,npx,1,npz+2)
                                                                   +L*sigma(0,npx,0,npz+3) ));
            sigma2(2,npx,1,npz)+= L*( sigma(2,npx,0,npz+1)+2.0*sigma(1,npx+1,1,npz)
                                      +L*(2.0*sigma(1,npx+1,0,npz+1)+sigma(0,npx+2,1,npz)
                                          +L*sigma(0,npx+2,0,npz+1) ));
            sigma2(1,npx,2,npz)+= L*( sigma(0,npx+1,2,npz)+2.0*sigma(1,npx,1,npz+1)
                                      +L*(2.0*sigma(0,npx+1,1,npz+1)+sigma(1,npx,0,npz+2)
                                          +L*sigma(0,npx+1,0,npz+2) ));
        }
    }
    // nx+nz=4
    sigma2(4,0,0,0)+= L*(4.0*sigma(3,1,0,0)
                         +L*(6.0*sigma(2,2,0,0)
                             +L*(4.0*sigma(1,3,0,0)
                                 +L*sigma(0,4,0,0) )));
    sigma2(0,0,4,0)+= L*(4.0*sigma(0,0,3,1)
                         +L*(6.0*sigma(0,0,2,2)
                             +L*(4.0*sigma(0,0,1,3)
                                 +L*sigma(0,0,0,4) )));
    sigma2(3,0,1,0)+= L*(sigma(3,0,0,1)+3.0*sigma(2,1,1,0)
                         +L*(3.0*sigma(2,1,0,1)+3.0*sigma(1,2,1,0)
                             +L*(3.0*sigma(1,2,0,1)+sigma(0,3,1,0)
                                 +L*sigma(0,3,0,1) )));
    sigma2(1,0,3,0)+= L*(sigma(0,1,3,0)+3.0*sigma(1,0,2,1)
                         +L*(3.0*sigma(0,1,2,1)+3.0*sigma(1,0,1,2)
                             +L*(3.0*sigma(0,1,1,2)+sigma(1,0,0,3)
                                 +L*sigma(0,1,0,3) )));
    sigma2(2,0,2,0)+= 2.0*L*(sigma(2,0,1,1)+sigma(1,1,2,0)
                         +L*(0.5*sigma(2,0,0,2)+2.0*sigma(1,1,1,1)+0.5*sigma(0,2,2,0)
                             +L*(sigma(1,1,0,2)+sigma(0,2,1,1)
                                 +0.5*L*sigma(0,2,0,2) )));
    sigma=sigma2;
}

void sc_Drift_Pass(STATvec & sigma, double L, double Ksc, unsigned Nint){

    const double ds = L/(double)Nint;

    // drift
    statVecMap_Drift(sigma,0.5*ds);
    for (unsigned i=0;i<Nint-1;i++){
        // SC kick
        statVecMap_SC_gaussian(sigma, ds*Ksc);
        // drift
        statVecMap_Drift(sigma,ds);
    }
    // SC kick
    statVecMap_SC_gaussian(sigma, ds*Ksc);
    // drift
    statVecMap_Drift(sigma,0.5*ds);
}

