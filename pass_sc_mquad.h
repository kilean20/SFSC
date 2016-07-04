#ifndef PASS_SC_MQUAD_H
#define PASS_SC_MQUAD_H

#include "statvec.h"

void statVecMap_mQuadKick(STATvec & sigma, double KL);
void sc_mQuadPass(STATvec & sigma, double  L, double K1, double Ksc, unsigned Nint);

#endif

