#include "line.h"
#include "element.h"
#include <vector>

using namespace std;

//=========================================================================
//
//                               Line  class
//
//=========================================================================
LINE::LINE ()
{
    Ncell = 0;
    Length = 0;
}

void LINE::Update ()
{
    double sPointer = 0.;
    for (short i = 0; i < Cell.size (); i++)
    {
        Cell[i].S = sPointer;
        sPointer = Cell[i].L + sPointer;
    }
    Ncell = Cell.size ();
    Length = sPointer;
}

void LINE::Append (ELEMENT elem)
{
    Cell.push_back (elem);
    elem.S=Length;
    Length += elem.L;
    Ncell+=1;
}
