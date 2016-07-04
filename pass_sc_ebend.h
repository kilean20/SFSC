#ifndef PASS_SC_EBEND_H
#define PASS_SC_EBEND_H

#include "statvec.h"

double dx_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom);
double dz_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom);
double dpx_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom);
double dpz_eBend(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double iRho, double  Ksc, double  denom);

void sc_eBendPass(STATvec & sigma, double  L, double Angle, double Ksc, unsigned Nint);

#endif

