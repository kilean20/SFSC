#include "global.h"
#include "pass_drift.h"

#include <cmath>
#include <vector>
using namespace std;

void DriftPass(vector<double> &x, double L)
{
    x[ps_]=sqrt( BETA2 *(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );
    const double lambda = L/x[ps_];
    x[x_] += x[px_] * lambda;
    x[z_] += x[pz_] * lambda;
    x[vt_] += (x[dE_]*BETA2 + 1.0) * lambda - L;
}


