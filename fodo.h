#include "elepass.h"

void fodo(Line & FODO)
{
    double ang = 2*datum::pi/12/2;
    double quadK = 5.0;

    double lBend = 0.4;
    double lQuad = 0.1;
    double lDrift = 0.2;

  Element * temp;
  for(int i=0;i<12;i++){
      temp = new DRIFT(0.5*lDrift);
      FODO.Append(temp);
      temp = new mBEND(lBend, ang);
      FODO.Append(temp);
      temp = new DRIFT(lDrift);
      FODO.Append(temp);
      temp = new mQUAD(lQuad, quadK);
      FODO.Append(temp);
      temp = new DRIFT(lDrift);
      FODO.Append(temp);
      temp = new mBEND(lBend, ang);
      FODO.Append(temp);
      temp = new DRIFT(lDrift);
      FODO.Append(temp);
      temp = new mQUAD(lQuad, -quadK);
      FODO.Append(temp);
      temp = new DRIFT(0.5*lDrift);
      FODO.Append(temp);
  }
}
