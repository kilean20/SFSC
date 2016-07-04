#include <iostream>
#include "statvec.h"
#include "pass_sc_drift.h"
#include "pass_sc_equad.h"

using namespace std;


int main()
{
    STATvec sigma;
    sigma(2,0,0,0) = 1e-4;
    sigma(0,2,0,0) = 1e-4;
    sigma(0,0,2,0) = 1e-5;
    sigma(0,0,0,2) = 1e-5;

    for(size_t dum=0;dum<100;dum++){
        for(size_t i=0;i<1;i++){
            //sc_Drift_Pass(sigma, 0.01, 0, 1);
            sc_eQuadPass(sigma, 0.5, 5.0, 0, 200);
            cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;
        }
        for(size_t i=0;i<1;i++){
            //sc_Drift_Pass(sigma, 0.01, 0, 1);
            sc_eQuadPass(sigma, 0.5, -5.0, 0, 200);
            cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;
        }
    }
    return 0;
}
