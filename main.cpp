#include <iostream>
#include <vector>
#include "efodo.h"
#include "line.h"
#include "global.h"
#include "statvec.h"



using namespace arma;


int main()
{
    cout.precision(15);
    //Define the ring
    LINE FODO;    efodo(FODO);    //FODO.Update();

    size_t nTurn=1E2;
    STATvec sigma;

    vector<STATvec> sigmaData();

        for(size_t i=0;i<nTurn;i++){
            for(size_t j=0;j<FODO.Ncell;j++) {FODO.Cell[j].sc_Pass(sigma);}
                cout << "x=" << sigma(1,0,0,0) << "  x2=" << sigma(2,0,0,0) << endl;
                sigmaData.push_back(sigma);
            }

    return 0;
}
