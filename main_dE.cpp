#include <iostream>
#include <vector>
#include "efodo.h"
#include "line.h"
#include "global.h"
#include <armadillo>
using namespace arma;

int main()
{
    cout.precision(15);
    //Define the ring
    Line FODO;    efodo(FODO);


    size_t nParticle=11, nTurn=1e7, nSkip=1e3;
    vector<double> x(10); fill(x.begin(), x.end(), 0.0);

    mat xData(nTurn/nSkip+1,2*nParticle);xData.zeros();
    mat zData(nTurn/nSkip+1,2*nParticle);zData.zeros();
    mat sData(nTurn/nSkip+1,2*nParticle);sData.zeros();
    mat spinData(nTurn/nSkip+1,2*nParticle);spinData.zeros();

    for(size_t n=0;n<nParticle;n++){
        fill(x.begin(), x.end(), 0.0);
        //initialization
        x[dE_]=((double)n-((double)nParticle-1.0)/2.0)*1.0e-4;
        for(size_t i=0;i<nTurn;i++){
            if (i%nSkip==0)
            {
                xData(i/nSkip,0+2*n)=x[x_];
                xData(i/nSkip,1+2*n)=x[px_];
                zData(i/nSkip,0+2*n)=x[z_];
                zData(i/nSkip,1+2*n)=x[pz_];
                sData(i/nSkip,0+2*n)=x[vt_];
                sData(i/nSkip,1+2*n)=x[dE_];
                spinData(i/nSkip,0+2*n)=x[Sx_];
                spinData(i/nSkip,1+2*n)=x[Ss_];
            }
            for(size_t j=0;j!=FODO.Ncell;j++) FODO.Cell[j].Pass(x);
        }
        xData(nTurn/nSkip,0+2*n)=x[x_];
        xData(nTurn/nSkip,1+2*n)=x[px_];
        zData(nTurn/nSkip,0+2*n)=x[z_];
        zData(nTurn/nSkip,1+2*n)=x[pz_];
        sData(nTurn/nSkip,0+2*n)=x[vt_];
        sData(nTurn/nSkip,1+2*n)=x[dE_];
        spinData(nTurn/nSkip,0+2*n)=x[Sx_];
        spinData(nTurn/nSkip,1+2*n)=x[Ss_];
    }
    xData.save("xDataTim",raw_ascii);
    zData.save("zDataTim",raw_ascii);
    sData.save("sDataTim",raw_ascii);
    spinData.save("spinDataTim",raw_ascii);
    return 0;
}
