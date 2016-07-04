#include "element.h"
#include "efodo.h"
#include "global.h"
#include <cmath>

void efodo(LINE & FODO)
{
    double ang = 2*M_PI/12/2;
    double quadK = 5.0;

    double lBend = 0.4;
    double lQuad = 0.1;
    double lDrift_ = 0.2;

  ELEMENT temp;
  for(int i=0;i<12;i++){
      temp.SetElem(DRIFT_,0.5*lDrift_);temp.Nint=30;
      FODO.Append(temp);
      temp.SetElem(eBEND_,lBend, ang);temp.Nint=120;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift_);temp.Nint=60;
      FODO.Append(temp);
      temp.SetElem(eQUAD_,lQuad, quadK);temp.Nint=30;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift_);temp.Nint=60;
      FODO.Append(temp);
      temp.SetElem(eBEND_,lBend, ang);temp.Nint=120;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift_);temp.Nint=60;
      FODO.Append(temp);
      temp.SetElem(eQUAD_,lQuad, -quadK);temp.Nint=30;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,0.5*lDrift_);temp.Nint=30;
      FODO.Append(temp);
  }
  const double cSpeed=299792458;
  const double revFreq = (BETA*cSpeed)/FODO.Length;

  temp.SetElem(RFcav_,1000.0, 3.0*revFreq, 0.0);
  FODO.Append(temp);
}
