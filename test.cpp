#include <iostream>
#include "statvec.h"
#include "sc_driftpass.h"
#include "sc_equadpass.h"
#include "sc_ebendpass.h"

using namespace std;


int main()
{
    STATvec sigma;
    sigma(2,0,0,0) = 1e-4;
    sigma(0,0,2,0) = 1e-5;
    //sigma(1,1,0,0) += sigma(2,0,0,0)+sigma(0,0,2,0);
    cout << sigma(1,1,0,0) << endl;

    sc_DriftPass(sigma, 0.5, 0.0, 4);
    cout << sigma(2,0,0,0) << endl;

    sc_eQuadPass(sigma, 0.5, 0.7, 0.0, 4);
    cout << sigma(2,0,0,0) << endl;

    sc_eBendPass(sigma, 0.4, 0.75, 0.0, 4);
    cout << sigma(2,0,0,0) << endl;



    return 0;
}
