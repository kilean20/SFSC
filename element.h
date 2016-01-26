#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include "statvec.h"

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
    short Type;

public:
    ELEMENT();
    ELEMENT(short Type, double l);
    ELEMENT(short Type, double l, double param);
    ELEMENT(short Type, double param1, double param2, double param3);

    void SetElem(short Type, double l);
    void SetElem(short Type, double l, double param);
    void SetElem(short Type, double param1, double param2, double param3);


    short Nint, Norder;
    bool FlagSpinTrack;
    double L, S, Ksc;
    void Pass (std::vector<double> &x);
    void sc_Pass (STATvec &sigma);
};


#endif
