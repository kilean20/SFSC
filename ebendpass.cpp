#include "ebendpass.h"

void eBendSpinPass2(vec &x, double L, double ANGLE, size_t Nint){
    x[s_]=0.0; // initialize s. s=0 is the element entrance
    const double iRHO=ANGLE/L;
    double hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;
    double dummy2;
    // initialize p_s : hard edge fringe field effect
    x[ps_]=hs*sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );

    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    //spin precession vector
    vec Lambda(3);
    double Cos, Sin;
    //-------------------------------
    //2nd Order Symplectic Integrator
    //-------------------------------
    for(size_t i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                            + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
        }
        if(i==0 || i==Nint-1)
        {
            //kick
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        //transverse drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //spin kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        Lambda[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
        dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
        Lambda[1] = -x[pz_]*dummy*iRHO/hs;
        Lambda[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
        //dummy=norm(Lambda);
        dummy=sqrt(Lambda[0]*Lambda[0]+Lambda[1]*Lambda[1]+Lambda[2]*Lambda[2]);
        dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
        if(dummy>1.0e-150){
            Lambda=Lambda/dummy;
            dummy*=lambda;
            Cos=cos(dummy);
            Sin=sin(dummy);
            dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
            x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
            x[Sx_]+=dummy;
        }
        else{
            dummy=2.0*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
            x[Ss_]+=2.0*lambda*(Lambda[2]*x[Sx_]-Lambda[0]*dummy2);
            x[Sx_]+=dummy;
        }
        //transverse drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        if(i>Nint-3)
        {
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        else
        {
            x[px_] += 2.0*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }

    //1st correction on extra or remaining legnth to the edge
    hs = (1.0 + x[x_]*iRHO );
    lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                    + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
    //kick
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    x[s_] += x[ps_]/(hs*hs) *lambda;
    x[vt_] += dummy *lambda;
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //spin kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    Lambda[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
    dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
    Lambda[1] = -x[pz_]*dummy*iRHO/hs;
    Lambda[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
    //dummy=norm(Lambda);
    dummy=sqrt(Lambda[0]*Lambda[0]+Lambda[1]*Lambda[1]+Lambda[2]*Lambda[2]);
    dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
    if(dummy>1.0e-150){
        Lambda=Lambda/dummy;
        dummy*=lambda;
        Cos=cos(dummy);
        Sin=sin(dummy);
        dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
        x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
        x[Sx_]+=dummy;
    }
    else{
        dummy=2.0*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
        x[Ss_]+=2.0*lambda*(Lambda[2]*x[Sx_]-Lambda[1]*dummy2);
        x[Sx_]+=dummy;
    }
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    x[s_] += x[ps_]/(hs*hs) *lambda;
    x[vt_] += dummy * lambda;
    //keep only deviation from the reference particle
    x[vt_] -= L;
}
void eBendSpinPass4(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendSpinPass6(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendSpinPass8(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendOrbitPass2(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendOrbitPass4(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendOrbitPass6(vec &x, double L, double ANGLE, size_t Ninit){}
void eBendOrbitPass8(vec &x, double L, double ANGLE, size_t Ninit){}

