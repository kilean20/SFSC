#ifndef ELEMENT_H
#define ELEMENT_H
#include <vector>
#include <cmath>
#include "global.h"
#include "statvec.h"

using namespace std;
//=========================================================================
//                            electric elements
//=========================================================================
//--------------------------------DRIFT------------------------------------
struct DRIFT
{
};
//--------------------------------eBEND------------------------------------
struct eBEND
{
    double Angle;
};
//--------------------------------eQUAD------------------------------------
struct eQUAD
{
    double K1;
};
//=========================================================================
//                             RF elements
//=========================================================================
//--------------------------------RFCAV------------------------------------
struct RFCAV
{
    double Vrf, Wrf, Phase;
};

//=========================================================================
//
//      Base class of elements - static polymorphism using union
//
//=========================================================================
class ELEMENT
{
private:
    union {
        DRIFT Drift;
        eBEND eBend;
        eQUAD eQuad;
        RFCAV RFcav;
    };
    size_t Type;

public:
    ELEMENT();
    ELEMENT(size_t Type, double l);
    ELEMENT(size_t Type, double l, double param);
    ELEMENT(size_t Type, double param1, double param2, double param3);

    void SetElem(size_t Type, double l);
    void SetElem(size_t Type, double l, double param);
    void SetElem(size_t Type, double param1, double param2, double param3);


    size_t Nint, Norder, FlagSpinTrack;
    double L, S, Ksc;
    void Pass (vector<double> &x);
    void sc_Pass (STATvec &sigma);
};


#endif
