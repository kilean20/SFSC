#ifndef PASS_MQUAD_H
#define PASS_MQUAD_H

#include <vector>

//mQuadSpinPass 2nd order method
//void mQuadSpinPass(std::vector<double> &x, double L, double K1, unsigned Nint);
//mQuadSpinPass High order composition method
//void mQuadSpinPass(std::vector<double> &x, double L, double K1, unsigned Nint, unsigned Norder);
//mQuadOrbitPass 2nd order method
void mQuadOrbitPass(std::vector<double> &x, double L, double K1, unsigned Nint);
//mQuadOrbitPass High order composition method
void mQuadOrbitPass(std::vector<double> &x, double L, double K1, unsigned Nint, unsigned Norder);

#endif

