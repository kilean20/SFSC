#include "element.h"
#include "ebendpass.h"
//=========================================================================
//
//                            Element classes
//
//=========================================================================
Element::Element()
{
    Norder=2;
    FlagSpinTrack=1;
}
//=========================================================================
//                           electric elements
//=========================================================================
//--------------------------------DRIFT------------------------------------
DRIFT::DRIFT(double l)
{
    L = l;
    TYPE=DRIFT_;
}
//-------------------------------DRIFT-Pass--------------------------------
void DRIFT::Pass (vec &x)
{
    x[ps_]=sqrt( BETA2 *(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );
    double lambda = L/x[ps_];
    x[x_] += x[px_] * lambda;
    x[z_] += x[pz_] * lambda;
    x[vt_] += (x[dE_]*BETA2 + 1.0) * lambda - L;
}
//---------------------------------eBEND-----------------------------------
eBEND::eBEND ( double l, double angle)
{
    TYPE=eBEND_;
    L=l; ANGLE=angle;
    Nint = ceil(ANGLE/0.02); //default Nint : fractional error is about (0.02)^(Norder+1)
}
//-------------------------------eBEND-Pass--------------------------------
void eBEND::Pass (vec &x)
{
    if(FlagSpinTrack){
        switch (Norder)
        {
        case 2:
            eBendSpinPass2(x,L,ANGLE,Nint);
            break;
        case 4:
            eBendSpinPass4(x,L,ANGLE,Nint);
            break;
        case 6:
            eBendSpinPass6(x,L,ANGLE,Nint);
            break;
        case 8:
            eBendSpinPass8(x,L,ANGLE,Nint);
            break;
        default:
            cout << "Err : Norder must 2,4,6 or 8" << endl;
        }
    }
    else{
        switch (Norder)
        {
        case 2:
            eBendOrbitPass2(x,L,ANGLE,Nint);
            break;
        case 4:
            eBendOrbitPass4(x,L,ANGLE,Nint);
            break;
        case 6:
            eBendOrbitPass6(x,L,ANGLE,Nint);
            break;
        case 8:
            eBendOrbitPass8(x,L,ANGLE,Nint);
            break;
        default:
            cout << "Err : Norder must 2,4,6 or 8" << endl;
        }
    }
}
//---------------------------------eQUAD-----------------------------------
eQUAD::eQUAD (double l, double k1)
{
    TYPE=eQUAD_;
    L=l; K1=k1;
    Nint = (L*L*abs(K1)/0.02); //default Nint :
}
//-------------------------------eQUAD-Pass--------------------------------
void eQUAD::Pass (vec &x)
{
    const double betaGammaK1 = BETAGAMMA*K1;
    double dummy = x[dE_]*BETA2 + 1.0 - 0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;
    double dummy2;
    //fringe field kick
    x[ps_]=sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );
    //timeStep
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // spin precession vector
    vec Lambda(3); //Lambda[0]=Lambda_x,  Lambda[1]=Lambda_s,  Lambda[2]=Lambda_z
    double Cos, Sin;

    //-------------------------------
    //2nd Order Symplectic Integrator
    //-------------------------------
    for(size_t i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }
        //drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //spin kick
        dummy = x[dE_]*BETA2 +1.0 -0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;
        dummy2 = BETA*(MDM+1.0/(1.0+GAMMA*dummy));
        Lambda[0] = betaGammaK1*( -0.5*x[x_]*EDM*dummy + x[ps_]*x[z_]*dummy2  );
        Lambda[1] = -betaGammaK1*( x[pz_]*x[x_] + x[px_]*x[z_]  )*dummy2;
        Lambda[2] = betaGammaK1*( 0.5*x[z_]*EDM*dummy + x[ps_]*x[x_]*dummy2  );
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
        //drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //kick
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        if(i==Nint-1)
        {
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }else
        {
            x[px_] -= dummy*K1*x[x_] *2.0*lambda;
            x[pz_] += dummy*K1*x[z_] *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }
    //fringe field kick
    x[vt_] -= L;
}
//=========================================================================
//                              RF elements
//=========================================================================
//--------------------------------RFCAV------------------------------------
//---------------------------RFCAV-constructor-----------------------------
RFCAV::RFCAV (double vRF, double fRF, double phase)
{
    VRF = vRF; // [Voltage]
    WRF = 2*datum::pi*fRF;
    PHASE = phase;
}
//-------------------------------RFCAV-Pass--------------------------------
void RFCAV::Pass (vec &x)
{
    double t = x[vt_]/(BETA*datum::c_0);
    double delta_dE = -VRF*1.0E-6/(BETA2*magicENERGY)*sin( WRF*t + PHASE);
    // sign corresponds to the negative charge
    x[dE_] += delta_dE;
    x[ps_]=sqrt( BETA2 *(x[dE_]+1.0/BETA2)*(x[dE_]+1.0/BETA2)-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );
}
