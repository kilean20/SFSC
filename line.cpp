#include "line.h"
//=========================================================================
//
//                               Line  class
//
//=========================================================================
Line::Line ()
{
    Ncell = 0;
    Length = 0;
}

void Line::Update ()
{
    double sPointer = 0.;
    for (size_t i = 0; i < Cell.size (); i++)
    {
        Cell[i]->S = sPointer;
        sPointer = Cell[i]->L + sPointer;
    }
    Ncell = Cell.size ();
    Length = sPointer;
}

void Line::Append (Element * elem)
{
    Cell.push_back (elem);
    Update ();
}
