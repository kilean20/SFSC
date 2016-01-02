#include <iostream>
#include <vector>
#include "efodo.h"
#include "line.h"

int main()
{
    cout.precision(15);
    //Define the ring
    Line FODO;    efodo(FODO);    FODO.Update();

    //add cavity
    double revFreq = (BETA*datum::c_0)/FODO.Length;
    cout << "length=" << FODO.Length << endl;
    cout << "revFreq=" << revFreq << endl;

    Element * temp;
    temp = new RFCAV(1000.0, 3.0*revFreq, 0.0);
    FODO.Append(temp);

    size_t nTurn=1e9;
    size_t nSkip=1e6;
    vec x(10); x.zeros();
    mat xData(nTurn/nSkip+1,2);xData.zeros();
    mat zData(nTurn/nSkip+1,2);zData.zeros();
    mat sData(nTurn/nSkip+1,2);sData.zeros();
    mat spinData(nTurn/nSkip+1,2);spinData.zeros();

    //x(dE_)=1.0e-4;
    for(size_t i=0;i<nTurn;i++){
        if (i%nSkip==0)
        {
            xData(i/nSkip,0)=x(x_);
            xData(i/nSkip,1)=x(px_);
            zData(i/nSkip,0)=x(z_);
            zData(i/nSkip,1)=x(pz_);
            sData(i/nSkip,0)=x(vt_);
            sData(i/nSkip,1)=x(dE_);
            spinData(i/nSkip,0)=x(Sx_);
            spinData(i/nSkip,1)=x(Ss_);
        }
        for(size_t j=0;j!=FODO.Ncell;j++) FODO.Cell[j]->Pass(x);
    }
    xData(nTurn/nSkip,0)=x(x_);
    xData(nTurn/nSkip,1)=x(px_);
    zData(nTurn/nSkip,0)=x(z_);
    zData(nTurn/nSkip,1)=x(pz_);
    sData(nTurn/nSkip,0)=x(vt_);
    sData(nTurn/nSkip,1)=x(dE_);
    spinData(nTurn/nSkip,0)=x(Sx_);
    spinData(nTurn/nSkip,1)=x(Ss_);

    xData.save("xData_dE0",raw_ascii);
    zData.save("zData_dE0",raw_ascii);
    sData.save("sData_dE0",raw_ascii);
    spinData.save("spinData_dE0",raw_ascii);
    return 0;
}
