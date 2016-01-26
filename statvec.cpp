#include "statvec.h"

//constructor
STATvec::STATvec():Zero(0){ for(short i=0;i<69;i++) StatVec[i]=0.0;}

// operator overloading (n)
double & STATvec::operator()(const short n)
{
    return StatVec[n];
}

// operator overloading (nx,npx,nz,npz)
double & STATvec::operator()(const short nx, const short npx, const short nz, const short npz)
{
    const short r=nx+npx+nz+npz;
    if(r > 4){ Zero=0.0; return Zero;}

    short nr = 0;
    double temp=1;
    for(short r_index=1;r_index<r;r_index++){
        temp=1;

        for(short j=1;j<4;j++){
            temp*=(r_index+j);
        }
        nr+=temp/6;
    }

    for(short r_index=r-nx+1;r_index<r+1;r_index++){
        temp=1;
        for(short j=1;j<3;j++){
            temp*=(r_index+j);
        }
        nr+=temp/2;
    }

    for(short r_index=r-nx-npx+1;r_index<r-nx+1;r_index++){
        temp=1;
        for(short j=1;j<2;j++){
            temp*=(r_index+j);
        }
        nr+=temp;
    }

    for(short r_index=r-nx-npx-nz+1;r_index<r-nx-npx+1;r_index++){
        nr+=1;
    }
    return StatVec[nr];
}

// operator overloading +
STATvec STATvec::operator+(STATvec rhs)
{
    STATvec out;
    for (short i=0;i<70;i++) out(i)=this->StatVec[i]+rhs(i);
    return out;
}

// operator overloading +=
void STATvec::operator+=(STATvec rhs)
{
    for (short i=0;i<70;i++) this->StatVec[i]+=rhs(i);
}

//inline
//STATvec operator+(const STATvec & lhs, const STATvec &rhs)
//{
//    STATvec lhsCopy=lhs,rhsCopy=rhs;
//    for (short i=0;i<70;i++) lhsCopy(i)+=rhsCopy(i);
//    return lhsCopy;
//}

//inline
//void operator+=(STATvec &lhs, const STATvec &rhs)
//{
//    STATvec rhsCopy=rhs;
//    for (short i=0;i<70;i++) lhs(i)+=rhsCopy(i);
//}
