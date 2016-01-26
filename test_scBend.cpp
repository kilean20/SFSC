#include <iostream>
#include "statvec.h"
#include "pass_sc_ebend.h"

using namespace std;


int main()
{
    STATvec sigma;
    sigma(2,0,0,0) = 1e-4;
    sigma(0,2,0,0) = 1e-4;
    sigma(0,0,2,0) = 1e-5;
    sigma(0,0,0,2) = 1e-5;


    sc_eBendPass(sigma, 2.0, 1.0, 0, 100000);
    cout << sigma(0,2,0,0) <<"\t"<< sigma(2,0,0,0) <<"\t"<< sigma(1,1,0,0) << endl;


    return 0;
}
