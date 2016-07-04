#ifndef LINE_H
#define LINE_H

#include <vector>
#include "element.h"

//=========================================================================
//
//                               Line  class
//
//=========================================================================
class LINE
{
public:
  LINE ();
  void Update ();
  void Append(ELEMENT elem);

  std::vector<ELEMENT> Cell;
  double Length;		//---length of reference orbit
  unsigned Ncell;			//---number of elements
};

#endif
