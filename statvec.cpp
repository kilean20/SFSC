#include "statvec.h"
#include <iostream>

//constructor
STATvec::STATvec():Zero(0){ for(unsigned i=0;i<69;i++) StatVec[i]=0.0;}

// operator overloading (n)
double & STATvec::operator()(const unsigned n)
{
    return StatVec[n];
}

// operator overloading (nx,npx,nz,npz)
double & STATvec::operator()(const unsigned nx, const unsigned npx, const unsigned nz, const unsigned npz)
{
    const unsigned r=nx+npx+nz+npz;
    if(r > 4){ Zero=0.0; return Zero;}

    unsigned nr = 0;
    double temp=1;
    for(unsigned r_index=1;r_index<r;r_index++){
        temp=1;

        for(unsigned j=1;j<4;j++){
            temp*=(r_index+j);
        }
        nr+=temp/6;
    }

    for(unsigned r_index=r-nx+1;r_index<r+1;r_index++){
        temp=1;
        for(unsigned j=1;j<3;j++){
            temp*=(r_index+j);
        }
        nr+=temp/2;
    }

    for(unsigned r_index=r-nx-npx+1;r_index<r-nx+1;r_index++){
        temp=1;
        for(unsigned j=1;j<2;j++){
            temp*=(r_index+j);
        }
        nr+=temp;
    }

    for(unsigned r_index=r-nx-npx-nz+1;r_index<r-nx-npx+1;r_index++){
        nr+=1;
    }
    return StatVec[nr];
}

// operator overloading +
STATvec STATvec::operator+(STATvec rhs)
{
    STATvec out;
    for (unsigned i=0;i<69;i++) out(i)=this->StatVec[i]+rhs(i);
    return out;
}

// operator overloading +=
void STATvec::operator+=(STATvec rhs)
{
    for (unsigned i=0;i<69;i++) this->StatVec[i]+=rhs(i);
}

//fill zeros
void STATvec::zeros()
{
    for(unsigned i=0;i<69;i++) StatVec[i]=0.0;
}
