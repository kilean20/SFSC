#include <iostream>
#include <vector>
#include "efodo.h"
#include "line.h"
#include "global.h"
#include "statvec.h"



using namespace std;


int main()
{
    cout.precision(15);
    //Define the ring
    LINE FODO;    efodo(FODO);    //FODO.Update();

    size_t nTurn=10;
    STATvec sigma;
    sigma(2,0,0,0)=1e-4;
    sigma(0,0,2,0)=1e-5;


    for(size_t i=0;i<nTurn;i++){
        for(size_t j=0;j<FODO.Ncell;j++) FODO.Cell[j].sc_Pass(sigma);
            cout << "x=" << sigma(1,0,0,0) << "  x2=" << sigma(2,0,0,0) << endl;
    }
    return 0;
}
