#include "global.h"
#include "pass_mquad.h"

#include <cmath>
#include <vector>
using namespace std;

//=============================================================================
//                       mQuadSpinPass 2nd order method
//=============================================================================
//=============================================================================
//                 mQuadSpinPass High order composition method
//=============================================================================
//=============================================================================
//                       mQuadOrbitPass 2nd order method
//=============================================================================
void mQuadOrbitPass(vector<double> &x, double L, double K1, unsigned Nint){

    //timeStep
    const double ds = 0.5*L/(double)Nint;
    double dummy;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= ds*K1*x[x_];
            x[pz_] += ds*K1*x[z_];
        }
        //drift
        dummy=sqrt(BETA2*(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);
        x[x_] += 2.0*ds*x[px_]/dummy;
        x[z_] += 2.0*ds*x[pz_]/dummy;
        x[vt_] += 2.0*ds*(x[dE_]*BETA2+1.0)/dummy;
        //kick
        if(i==Nint-1)
        {
            x[px_] -= ds*K1*x[x_];
            x[pz_] += ds*K1*x[z_];
        }else
        {
            x[px_] -= 2.0*ds*K1*x[x_];
            x[pz_] += 2.0*ds*K1*x[z_];
        }
    }
}
//=============================================================================
//                eQuadOrbitPass High order composition method
//=============================================================================
void mQuadOrbitPass(vector<double> &x, double L, double K1, unsigned Nint, unsigned Norder){
    // initialize step size
    vector<double> R;
    unsigned nR;
    switch (Norder) {
    case 2:
        R.push_back(1.0);
        nR = 1;
        break;
    case 4:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    case 6:
        for(unsigned i=0;i<7;i++)
        R.push_back(R6[i]);
        nR = 7;
        break;
    default:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    }

    //timeStep
    const double ds = 0.5*L/(double)Nint;
    double dummy;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= R[0]*ds*K1*x[x_];
            x[pz_] += R[0]*ds*K1*x[z_];
        }
        //drift
        dummy=sqrt(BETA2*(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);
        x[x_] += 2.0*R[0]*ds*x[px_]/dummy;
        x[z_] += 2.0*R[0]*ds*x[pz_]/dummy;
        x[vt_] += 2.0*R[0]*ds*(x[dE_]*BETA2+1.0)/dummy;
        for(unsigned r=1; r<nR; r++)
        {
            //kick
            x[px_] -= (R[r]+R[r-1])*ds*K1*x[x_];
            x[pz_] += (R[r]+R[r-1])*ds*K1*x[z_];
            //drift
            dummy=sqrt(BETA2*(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);
            x[x_] += 2.0*R[r]*ds*x[px_]/dummy;
            x[z_] += 2.0*R[r]*ds*x[pz_]/dummy;
            x[vt_] += 2.0*R[r]*ds*(x[dE_]*BETA2+1.0)/dummy;
        }
        //kick
        if(i==Nint-1)
        {
            x[px_] -= R[nR-1]*ds*K1*x[x_];
            x[pz_] += R[nR-1]*ds*K1*x[z_];
        }else
        {
            x[px_] -= 2.0*R[nR-1]*ds*K1*x[x_];
            x[pz_] += 2.0*R[nR-1]*ds*K1*x[z_];
        }
    }
}

