#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include "statvec.h"


//=========================================================================
//
//                               Elements
//
//=========================================================================
//--------------------------------DRIFT------------------------------------
struct DRIFT
{
};
//=========================================================================
//                            electric elements
//=========================================================================
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
//                            magnetic elements
//=========================================================================
//--------------------------------mBEND------------------------------------
struct mBEND
{
    double Angle;
};
//--------------------------------mQUAD------------------------------------
struct mQUAD
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
        DRIFT Drift_;
        eBEND eBend;
        eQUAD eQuad;
        mBEND mBend;
        mQUAD mQuad;
        RFCAV RFcav;
    };
    unsigned Type;

public:
    ELEMENT();
    ELEMENT(unsigned Type, double l);
    ELEMENT(unsigned Type, double l, double param);
    ELEMENT(unsigned Type, double param1, double param2, double param3);

    void SetElem(unsigned Type, double l);
    void SetElem(unsigned Type, double l, double param);
    void SetElem(unsigned Type, double param1, double param2, double param3);


    unsigned Nint, Norder;
    bool FlagSpinTrack;
    double L, S, Ksc;
    void Pass (std::vector<double> &x);
    void sc_Pass (STATvec &sigma);
};


#endif
