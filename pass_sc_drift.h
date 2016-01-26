#ifndef PASS_SC_DRIFT_H
#define PASS_SC_DRIFT_H

#include "statvec.h"

double dx_Drift(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dz_Drift(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dpx_Drift(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dpz_Drift(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double  Ksc, double  denom);
void sc_DriftPass(STATvec & sigma, double  L, double  Ksc, short Nint);

#endif // DRIFTPASS_H

