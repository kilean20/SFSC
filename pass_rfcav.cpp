#include "global.h"

#include <vector>
#include <cmath>

using namespace std;

void RFcavPass(vector<double> &x, double VRF, double WRF, double PHASE)
{
    const double cSpeed=299792458;
    const double t = x[vt_]/(BETA*cSpeed);
    double delta_dE = -VRF*1.0E-6/(BETA2*magicENERGY)*sin( WRF*t + PHASE);
    // sign corresponds to the negative charge
    x[dE_] += delta_dE;
}

