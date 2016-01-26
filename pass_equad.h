#ifndef PASS_EQUAD_H
#define PASS_EQUAD_H

#include <vector>

//eQuadSpinPass 2nd order method
void eQuadSpinPass(std::vector<double> &x, double L, double K1, short Nint);
//eQuadSpinPass High order composition method
void eQuadSpinPass(std::vector<double> &x, double L, double K1, short Nint, short Norder);
//eQuadOrbitPass 2nd order method
void eQuadOrbitPass(std::vector<double> &x, double L, double K1, short Nint);
//eQuadOrbitPass High order composition method
void eQuadOrbitPass(std::vector<double> &x, double L, double K1, short Nint, short Norder);

#endif  //EQUADPASS

