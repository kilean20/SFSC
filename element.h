#ifndef ELEMENT_H
#define ELEMENT_H

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <include/armadillo>
#else
    #include <armadillo>
#endif

#include <vector>
#include <cmath>
#include "global.h"

using namespace std;
using namespace arma;

//=========================================================================
//
//                            Element classes
//
//=========================================================================
class Element
{
public:
  size_t Nint, TYPE, Norder, FlagSpinTrack;
  double L, S;
  virtual void Pass (vec &x)=0;
  Element();
};
//=========================================================================
//                           electric elements
//=========================================================================
//--------------------------------DRIFT------------------------------------
class DRIFT:public Element
{
public:
  DRIFT (double length);
  void Pass(vec &x);
};
//--------------------------------eBEND------------------------------------
class eBEND:public Element
{
public:
  eBEND (double l, double angle);
  void Pass(vec &x);
  double ANGLE;
};
//--------------------------------eQUAD------------------------------------
class eQUAD:public Element
{
public:
  eQUAD (double l, double k1);
  void Pass(vec &x);
  double K1;
};
//=========================================================================
//                           magnetic elements
//=========================================================================
////--------------------------------mBEND------------------------------------
//class mBEND:public Element
//{
//public:
//  mBEND (double l, double angle);
//  void Pass (vec &x);
//  double ANGLE, RHO;
//};
////------------------------------mBEND_EDGE---------------------------------
//class mBEND_EDGE:public Element
//{
//public:
//  mBEND_EDGE (double edge, double rho, double fint, double gap, bool inout);
//  mBEND_EDGE (double edge, double rho);
//  void Pass (vec &x);
//  double EDGE, RHO, FINT, GAP;
//  bool InOut;
//};
////---------------------------------mQUAD-----------------------------------
//class mQUAD:public Element
//{
//public:
// mQUAD (double l, double k1);
// void Pass (vec &x);
// double K1;
//};
//=========================================================================
//                              RF elements
//=========================================================================
//--------------------------------RFCAV------------------------------------
class RFCAV:public Element
{
public:
  RFCAV (double vRF, double wRF, double phase);
  void Pass (vec &x);
  double VRF, WRF, PHASE;
};

#endif
