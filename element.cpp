#include "element.h"
#include "global.h"
#include "pass_drift.h"
#include "pass_ebend.h"
#include "pass_equad.h"
#include "pass_rfcav.h"
#include "pass_sc_drift.h"
#include "pass_sc_ebend.h"
#include "pass_sc_mquad.h"
//#include "pass_sc_mbend.h"

#include <vector>
#include <cmath>
using namespace std;

//=========================================================================
//
//                            Element classes
//
//=========================================================================
//-------------------------------constructor0------------------------------
ELEMENT::ELEMENT():Type(0),L(0),Ksc(0){Norder=2; FlagSpinTrack=1;}
//-------------------------------constructor1------------------------------
ELEMENT::ELEMENT(unsigned type, double l):Type(type),L(l),Ksc(0){
    Nint = 4;
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case DRIFT_:
        DRIFT d;
        Drift_=d;
        break;
    default:
        break;
    }
}
//------------------------------constructor2--------------------------------
ELEMENT::ELEMENT(unsigned type, double l,double param):Type(type),L(l),Ksc(0){
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case eBEND_:
        eBEND b;  b.Angle=param;
        eBend=b;
        Nint = ceil(param/0.02); //default Nint
        break;
    case eQUAD_:
        eQUAD q;  q.K1=param;
        eQuad=q;
        Nint = (L*L*abs(param)/0.02); //default Nint
        break;
    default:
        break;
    }
}
//------------------------------constructor3--------------------------------
ELEMENT::ELEMENT(unsigned type, double v,double f,double phi):Type(type),L(0),Ksc(0){
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case RFcav_:
        RFCAV rf;   rf.Vrf=v;   rf.Wrf=2*M_PI*f; rf.Phase=phi;
        RFcav=rf;
    default:
        break;
    }
}
//=================================SetElem=================================
//---------------------------------SetElem1--------------------------------
void
ELEMENT::SetElem(unsigned type, double l){
    Type=type; L=l;
    switch (type){
    case DRIFT_:
        DRIFT d;
        Drift_=d;
        break;
    default:
        break;
    }
}
//---------------------------------SetElem2--------------------------------
void
ELEMENT::SetElem(unsigned type, double l, double param){
    Type=type; L=l;
    switch (type){
    case eBEND_:
        eBEND ebend;  ebend.Angle=param;
        eBend=ebend;
        Nint = ceil(param/0.02); //default Nint
        break;
    case eQUAD_:
        eQUAD equad;  equad.K1=param;
        eQuad=equad;
        Nint = (L*L*abs(param)/0.02); //default Nint
        break;
    case mBEND_:
        mBEND mbend;  mbend.Angle=param;
        mBend=mbend;
        Nint = ceil(param/0.02); //default Nint
        break;
    case mQUAD_:
        mQUAD mquad;  mquad.K1=param;
        mQuad=mquad;
        Nint = (L*L*abs(param)/0.02); //default Nint
        break;
    default:
        break;
    }
}
//---------------------------------SetElem3--------------------------------
void
ELEMENT::SetElem(unsigned type, double v,double f,double phi){
    Type=type;
    switch (type){
    case RFcav_:
        RFCAV rf;   rf.Vrf=v;   rf.Wrf=2*M_PI*f; rf.Phase=phi;
        RFcav=rf;
    default:
        break;
    }
}
//=========================================================================
//                              Pass method
//=========================================================================
void
ELEMENT::Pass (vector<double> &x)
{
    switch (Type) {
    case DRIFT_:
        Drift_Pass(x, L);
        break;
    case eBEND_:
        if (FlagSpinTrack)
            if (Norder==2)
                eBendSpinPass(x, L, eBend.Angle, Nint);
            else
                eBendSpinPass(x, L, eBend.Angle, Nint, Norder);
        else
            if (Norder==2)
                eBendOrbitPass(x, L, eBend.Angle, Nint);
            else
                eBendOrbitPass(x, L, eBend.Angle, Nint, Norder);
        break;
    case eQUAD_:
        if (FlagSpinTrack)
            if (Norder==2)
                eQuadSpinPass(x, L, eQuad.K1, Nint);
            else
                eQuadSpinPass(x, L, eQuad.K1, Nint, Norder);
        else
            if (Norder==2)
                eQuadOrbitPass(x, L, eQuad.K1, Nint);
            else
                eQuadOrbitPass(x, L, eQuad.K1, Nint, Norder);
        break;
    case mQUAD_:
//        if (FlagSpinTrack)
//            if (Norder==2)
//                mQuadSpinPass(x, L, eQuad.K1, Nint);
//            else
//                mQuadSpinPass(x, L, eQuad.K1, Nint, Norder);
//        else
            if (Norder==2)
                mQuadOrbitPass(x, L, eQuad.K1, Nint);
            else
                mQuadOrbitPass(x, L, eQuad.K1, Nint, Norder);
        break;
    case RFcav_:
        RFcavPass(x, RFcav.Vrf, RFcav.Wrf, RFcav.Phase);
        break;
    default:
        break;
    }
}
//=========================================================================
//             Space Charge Pass Method of Statistical Moments
//=========================================================================
void
ELEMENT::sc_Pass (STATvec &sigma)
{
    switch (Type) {
    case DRIFT_:
        sc_Drift_Pass(sigma, L, Ksc, Nint);
        break;
    case eBEND_:
        sc_eBendPass(sigma, L, eBend.Angle, Ksc, Nint);
        break;
    case eQUAD_:
        sc_eQuadPass(sigma, L, eQuad.K1, Ksc, Nint);
        break;
    case mBEND_:
        //sc_mBendPass(sigma, L, mBend.Angle, Ksc, Nint);
        break;
    case mQUAD_:
        sc_mQuadPass(sigma, L, mQuad.K1, Ksc, Nint);
        break;
    case RFcav_:
        break;
    default:
        break;
    }
}
