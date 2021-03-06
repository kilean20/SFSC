#ifndef PASS_EBEND_H
#define PASS_EBEND_H

#include <vector>


//eBendSpinPass 2nd order method
void eBendSpinPass(std::vector<double> &x, double L, double ANGLE, unsigned Nint);
//eBendSpinPass High order composition method
void eBendSpinPass(std::vector<double> &x, double L, double ANGLE, unsigned Nint, unsigned Norder);
//eBendOrbitPass 2nd order method
void eBendOrbitPass(std::vector<double> &x, double L, double ANGLE, unsigned Nint);
//eBendOrbitPass High order composition method
void eBendOrbitPass(std::vector<double> &x, double L, double ANGLE, unsigned Nint, unsigned Norder);

#endif // EBENDPASS_H

