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

    double & operator()(const short n);
    double & operator()(const short nx, const short npx, const short nz, const short npz);

    STATvec operator+(STATvec rhs);
    void operator+=(STATvec rhs);

};

#endif
