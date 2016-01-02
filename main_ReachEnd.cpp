#include <iostream>
#include "element.h"
#include "line.h"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <include/armadillo>
#else
    #include <armadillo>
#endif

using namespace std;

int main()
{
    cout.precision(20);
    eBEND testBend(0.4,datum::pi/2);

    //test
    vec x(10);

    size_t nParticle = 1000;
    vec ds(nParticle);
    ds.zeros();
    for (size_t n=0;n<nParticle;n++)
    {
        x.zeros();
        x(x_)=randn()*1e-3;
        x(z_)=randn()*1e-3;
        x(px_)=randn()*1e-3;
        x(pz_)=randn()*1e-3;
        x(dE_)=randn()*1e-3;
        x(ps_)=sqrt( BETA2 *(x(dE_)+1.0/BETA2)*(x(dE_)+1.0/BETA2)-1.0/BETAGAMMA2 - x(px_)*x(px_) - x(pz_)*x(pz_) );
        testBend.Pass(x);
        ds(n)=x(s_)-0.4;
    }
    cout << stddev(ds)  << endl;


    x.zeros();
    x(ps_)=1;
    testBend.Pass(x);
    cout << x(s_)-0.4 << endl;

    x.zeros();
    x(x_)=1e-3;
    x(dE_)=1e-3;
    x(ps_)=sqrt( BETA2 *(x(dE_)+1.0/BETA2)*(x(dE_)+1.0/BETA2)-1.0/BETAGAMMA2 - x(px_)*x(px_) - x(pz_)*x(pz_) );
    testBend.Pass(x);
    cout << x(s_)-0.4 << endl;


    x.zeros();
    x(x_)=2e-3;
    x(z_)=2e-3;
    x(px_)=2e-3;
    x(pz_)=2e-3;
    x(dE_)=2e-3;
    x(ps_)=sqrt( BETA2 *(x(dE_)+1.0/BETA2)*(x(dE_)+1.0/BETA2)-1.0/BETAGAMMA2 - x(px_)*x(px_) - x(pz_)*x(pz_) );
    testBend.Pass(x);
    cout << x(s_)-0.4 << endl;
	return 0;
}

