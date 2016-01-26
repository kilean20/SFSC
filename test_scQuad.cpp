#include <iostream>
#include "statvec.h"
#include "pass_sc_equad.h"
#include "pass_sc_drift.h"

using namespace std;


int main()
{
    STATvec sigma;
    sigma(2,0,0,0) = 1e-4;
    sigma(0,2,0,0) = 1e-4;
    sigma(0,0,2,0) = 1e-5;
    sigma(0,0,0,2) = 1e-5;


    //sc_eQuadPass(sigma, 0.75, 0.0, 0.0, 10);
    sc_eQuadPass(sigma, 0.75, 0.0, 10);
    //sc_DriftPass(sigma, 0.75, 0, 10);
    cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;


    return 0;
}
