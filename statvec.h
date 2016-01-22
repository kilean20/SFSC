#ifndef STATVEC_H
#define STATVEC_H
#include <vector>

using namespace std;
//=========================================================================
//                    Statistical Moment container
//=========================================================================
class STATvec
{
private:
    double StatVec[69];
    double zero;
public:
    STATvec(){ for(size_t i=0;i<69;i++) StatVec[i]=0.0; zero=0.0;}

    //double operator()(const size_t nx, const size_t npx, const size_t nz, const size_t npz);
    double & operator()(const size_t nx, const size_t npx, const size_t nz, const size_t npz);
};

//inline double STATvec::operator()(const size_t nx, const size_t npx, const size_t nz, const size_t npz)
inline double & STATvec::operator()(const size_t nx, const size_t npx, const size_t nz, const size_t npz)
{
    const size_t r=nx+npx+nz+npz;
    if(r > 4){ zero=0.0; return zero;}

    size_t nr = 0;
    double temp=1;
    for(size_t r_index=1;r_index<r;r_index++){
        temp=1;

        for(size_t j=1;j<4;j++){
            temp*=(r_index+j);
        }
        nr+=temp/6;
    }

    for(size_t r_index=r-nx+1;r_index<r+1;r_index++){
        temp=1;
        for(size_t j=1;j<3;j++){
            temp*=(r_index+j);
        }
        nr+=temp/2;
    }

    for(size_t r_index=r-nx-npx+1;r_index<r-nx+1;r_index++){
        temp=1;
        for(size_t j=1;j<2;j++){
            temp*=(r_index+j);
        }
        nr+=temp;
    }

    for(size_t r_index=r-nx-npx-nz+1;r_index<r-nx-npx+1;r_index++){
        nr+=1;
    }

    return StatVec[nr];
}
#endif
