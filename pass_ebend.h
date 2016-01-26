#ifndef PASS_EBEND_H
#define PASS_EBEND_H

#include <vector>


//eBendSpinPass 2nd order method
void eBendSpinPass(std::vector<double> &x, double L, double ANGLE, short Nint);
//eBendSpinPass High order composition method
void eBendSpinPass(std::vector<double> &x, double L, double ANGLE, short Nint, short Norder);
//eBendOrbitPass 2nd order method
void eBendOrbitPass(std::vector<double> &x, double L, double ANGLE, short Nint);
//eBendOrbitPass High order composition method
void eBendOrbitPass(std::vector<double> &x, double L, double ANGLE, short Nint, short Norder);

#endif // EBENDPASS_H

