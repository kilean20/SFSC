#include "element.h"
#include "efodoLinac.h"
#include "global.h"

void efodoLinac(LINE & FODO)
{
    double quadK = 7.64;

    double lQuad = 0.25;
    double lDrift = 0.75;

    ELEMENT temp;
    temp.SetElem(DRIFT_,0.5*lDrift);temp.Nint=250;
    FODO.Append(temp);
    temp.SetElem(eQUAD_,lQuad, quadK);temp.Nint=300;
    FODO.Append(temp);
    temp.SetElem(DRIFT_,lDrift);temp.Nint=500;
    FODO.Append(temp);
    temp.SetElem(eQUAD_,lQuad, -quadK);temp.Nint=300;
    FODO.Append(temp);
    temp.SetElem(DRIFT_,0.5*lDrift);temp.Nint=250;
    FODO.Append(temp);
}
