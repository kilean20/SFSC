#ifndef PASS_SC_DRIFT_H
#define PASS_SC_DRIFT_H

#include "statvec.h"

void statVecMap_Drift(STATvec & sigma, double L);
void sc_Drift_Pass(STATvec & sigma, double L, double Ksc, unsigned Nint);

#endif

