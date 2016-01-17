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
    LINE FODO;    efodo(FODO);

    size_t nTurn=1e9;
    size_t nSkip=1e6;
    vector<double> x(10); std::fill(x.begin(), x.end(), 0.0);

    mat xData(nTurn/nSkip+1,2);xData.zeros();
    mat zData(nTurn/nSkip+1,2);zData.zeros();
    mat sData(nTurn/nSkip+1,2);sData.zeros();
    mat spinData(nTurn/nSkip+1,2);spinData.zeros();

    x[dE_]=-5.0e-4;
    for(size_t i=0;i<nTurn;i++){
        if (i%nSkip==0)
        {
            xData(i/nSkip,0)=x[x_];
            xData(i/nSkip,1)=x[px_];
            zData(i/nSkip,0)=x[z_];
            zData(i/nSkip,1)=x[pz_];
            sData(i/nSkip,0)=x[vt_];
            sData(i/nSkip,1)=x[dE_];
            spinData(i/nSkip,0)=x[Sx_];
            spinData(i/nSkip,1)=x[Ss_];
        }
        for(size_t j=0;j!=FODO.Ncell;j++) FODO.Cell[j].Pass(x);
    }
    xData(nTurn/nSkip,0)=x[x_];
    xData(nTurn/nSkip,1)=x[px_];
    zData(nTurn/nSkip,0)=x[z_];
    zData(nTurn/nSkip,1)=x[pz_];
    sData(nTurn/nSkip,0)=x[vt_];
    sData(nTurn/nSkip,1)=x[dE_];
    spinData(nTurn/nSkip,0)=x[Sx_];
    spinData(nTurn/nSkip,1)=x[Ss_];

    xData.save("xData_dE_5",raw_ascii);
    zData.save("zData_dE_5",raw_ascii);
    sData.save("sData_dE_5",raw_ascii);
    spinData.save("spinData_dE_5",raw_ascii);
    return 0;
}
