#ifndef PASS_SC_EQUAD_H
#define PASS_SC_EQUAD_H

#include "statvec.h"

double dx_eQuad(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dz_eQuad(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dpx_eQuad(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);
double dpz_eQuad(short nx, short npx, short nz, short npz, STATvec & sigma, double  r, double K1, double  Ksc, double  denom);

void sc_eQuadPass(STATvec & sigma, double  L, double K1, double Ksc, short Nint);
void sc_eQuadPass(STATvec & sigma, double  L, double Ksc, short Nint);

#endif

