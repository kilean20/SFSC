#ifndef PASS_SC_EQUAD_H
#define PASS_SC_EQUAD_H

#include "statvec.h"

double dx_eQuad(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dz_eQuad(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dpx_eQuad(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dpz_eQuad(unsigned nx, unsigned npx, unsigned nz, unsigned npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);

void sc_eQuadPass(STATvec & sigma, double  L, double K1, double Ksc, unsigned Nint);

#endif

