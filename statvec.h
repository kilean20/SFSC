#ifndef STATVEC_H
#define STATVEC_H

//=========================================================================
//                    Statistical Moment container
//=========================================================================
class STATvec
{
     double Zero;
     double StatVec[69];
public:
    STATvec();

    double & operator()(const unsigned n);
    double & operator()(const unsigned nx, const unsigned npx, const unsigned nz, const unsigned npz);

    STATvec operator+(STATvec rhs);
    void operator+=(STATvec rhs);
    void zeros();

};

#endif
