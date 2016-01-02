#ifndef EBENDPASS_H
#define EBENDPASS_H


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <include/armadillo>
#else
    #include <armadillo>
#endif
#include "global.h"
using namespace arma;
using namespace std;


void eBendSpinPass2(vec &x, double L, double ANGLE, size_t Ninit);
void eBendSpinPass4(vec &x, double L, double ANGLE, size_t Ninit);
void eBendSpinPass6(vec &x, double L, double ANGLE, size_t Ninit);
void eBendSpinPass8(vec &x, double L, double ANGLE, size_t Ninit);
void eBendOrbitPass2(vec &x, double L, double ANGLE, size_t Ninit);
void eBendOrbitPass4(vec &x, double L, double ANGLE, size_t Ninit);
void eBendOrbitPass6(vec &x, double L, double ANGLE, size_t Ninit);
void eBendOrbitPass8(vec &x, double L, double ANGLE, size_t Ninit);

#endif // EBENDPASS_H

