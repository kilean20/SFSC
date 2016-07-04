#ifndef PASS_RK4_SC_DRIFT_H
#define PASS_RK4_SC_DRIFT_H

#include "statvec.h"

double dx_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dz_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dpx_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom);
double dpz_Drift_RK4(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double  Ksc, double  denom);
void sc_Drift_Pass_RK4(STATvec & sigma, double  L, double  Ksc, unsigned Nint);

#endif // DRIFTPASS_H

