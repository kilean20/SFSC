#include <iostream>
#include "statvec.h"
#include "sc_driftpass.h"

using namespace std;


int main()
{
    STATvec sigma;
    // operator overloading "(nx,npx,nz,npz)"
    sigma(2,0,0,0) = 1e-4;
    //sigma(0,2,0,0) = 1e-4;
    sigma(0,0,2,0) = 1e-5;
    sigma(0,0,0,1) = 5e-8;
    // operator overloading "(n)"
    sigma(2) = 3e-8;
    cout << sigma(0) <<"\t"<< sigma(1) <<"\t"<< sigma(2) << endl;
    cout << sigma(0,0,0,1) <<"\t"<< sigma(0,0,1,0) <<"\t"<< sigma(0,1,0,0) << endl;

    // operator overloading "+"
    STATvec temp=sigma+sigma;
    cout << sigma(0) <<"\t"<< sigma(1) <<"\t"<< sigma(2) << endl;
    cout << temp(0) <<"\t"<< temp(1) <<"\t"<< temp(2) << endl;

    // operator overloading "+="
    sigma+=sigma+sigma;
    cout << sigma(0) <<"\t"<< sigma(1) <<"\t"<< sigma(2) << endl;

    // reset initial
    sigma(0,0,0,1) = 0.0;
    sigma(2) = 0.0;
    cout << sigma(2,0,0,0) <<"\t"<< sigma(0,0,2,0) <<"\t"<< sigma(0,1,0,0) << endl;


    sc_DriftPass(sigma, 0.5, 1e-8, 20);
    cout << sigma(2,0,0,0) <<"\t"<< sigma(0,0,2,0) <<"\t"<< sigma(0,1,0,0) << endl;

    sc_eQuadPass(sigma, 0.5, 0.7, 0.0, 4);
    cout << sigma(2,0,0,0) <<"\t"<< sigma(0,0,2,0) <<"\t"<< sigma(0,1,0,0) << endl;

    sc_eBendPass(sigma, 0.4, 0.75, 0.0, 4);
    cout << sigma(2,0,0,0) <<"\t"<< sigma(0,0,2,0) <<"\t"<< sigma(0,1,0,0) << endl;



    return 0;
}
