#ifndef LINE_H
#define LINE_H

#include <vector>
#include "global.h"
#include "element.h"

using namespace std;
//=========================================================================
//
//                               Line  class
//
//=========================================================================
class Line
{
public:
  Line ();
  void Update ();
  void Append(Element * x);

  vector<Element * > Cell;
  double Length;		//---length of reference orbit
  size_t Ncell;			//---number of elements
};

#endif
