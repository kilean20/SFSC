#include "elepass.h"
//=========================================================================
//
//                            Element classes
//
//=========================================================================
Element::Element()
{
    //----def optim time step coeffi
    R[0]= 0.74167036435061295344822780;
    R[1]=-0.40910082580003159399730010;
    R[2]= 0.19075471029623837995387626;
    R[3]=-0.57386247111608226665638773;
    R[4]= 0.29906418130365592384446354;
    R[5]= 0.33462491824529818378495798;
    R[6]= 0.31529309239676659663205666;
    R[7]=-0.79688793935291635401978884;
    for(size_t r=1;r<8;r++)  R[7+r]=R[7-r];
}

//=========================================================================
//                           electric elements
//=========================================================================
//--------------------------------DRIFT------------------------------------
DRIFT::DRIFT( double l)
{
    TYPE=0;
    L = l;
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
    TYPE=1;
    L=l; ANGLE=angle; RHO=l/angle;
    //Nint = ceil(0.1*L/RHO);  // default number of integration step
    size_t Nkick = 300;
    Nint = ceil( (float)Nkick/(float)15);
}
//-------------------------------eBEND-Pass--------------------------------
void eBEND::Pass (vec &x)
{
    x[s_]=0.0; // initialize s. s=0 is the element entrance
    double hs, iRHO=1.0/RHO;
    hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;
    double dummy2;
    //fringe field kick
    x[ps_]=hs*sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );
    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;
    //spin precession vector
    vec Lambda(3);
    double Cos, Sin;
    //-------------------------------
    //8th Order Symplectic Integrator
    //-------------------------------
    for(size_t i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                            + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
        }
        //kick
        x[px_] += R[0]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
        x[s_] += x[ps_]/(hs*hs) *R[0]*lambda;
        x[vt_] += dummy *R[0]*lambda;
        //transverse drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
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
            dummy*=R[0]*lambda;
            Cos=cos(dummy);
            Sin=sin(dummy);
            dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
            x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
            x[Sx_]+=dummy;
        }
        else{
            dummy=2.0*R[0]*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
            x[Ss_]+=2.0*R[0]*lambda*(Lambda[2]*x[Sx_]-Lambda[0]*dummy2);
            x[Sx_]+=dummy;
        }
        //transverse drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;

        for(size_t r=1; r<15; r++)
        {
            //kick
            hs = (1.0 + x[x_]*iRHO );
            dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
            x[px_] += (R[r]+R[r-1])*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy * (R[r]+R[r-1])*lambda;
            //transverse drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
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
                dummy*=R[r]*lambda;
                Cos=cos(dummy);
                Sin=sin(dummy);
                dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
                x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
                x[Sx_]+=dummy;
            }
            else{
                dummy=2.0*R[r]*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
                x[Ss_]+=2.0*R[r]*lambda*(Lambda[2]*x[Sx_]-Lambda[0]*dummy2);
                x[Sx_]+=dummy;
            }
            //transverse drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
        }
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        x[px_] += R[14]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
        x[s_] += x[ps_]/(hs*hs) *R[14]*lambda;
        x[vt_] += dummy *R[14]*lambda;
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
    x[vt_] += dummy * lambda - L;

    // vt deviation
    //    x[vt_] -= L;
}
//---------------------------------eQUAD-----------------------------------
eQUAD::eQUAD (double l, double k1)
{
    TYPE=2;
    L=l; K1=k1;
    //Nint=ceil(0.01*L*L*K1);
    size_t Nkick = 75;
    Nint = ceil( (float)Nkick/(float)15);
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
    //8th Order Symplectic Integrator
    //-------------------------------
    for(size_t i=0; i<Nint ;i++)
    {
        //kick
        x[px_] -= dummy*K1*x[x_] *R[0]*lambda;
        x[pz_] += dummy*K1*x[z_] *R[0]*lambda;
        x[vt_] += dummy *R[0]*lambda;
        //drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
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
            dummy*=R[0]*lambda;
            Cos=cos(dummy);
            Sin=sin(dummy);
            dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
            x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
            x[Sx_]+=dummy;
        }
        else{
            dummy=2.0*R[0]*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
            x[Ss_]+=2.0*R[0]*lambda*(Lambda[2]*x[Sx_]-Lambda[0]*dummy2);
            x[Sx_]+=dummy;
        }
        //drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
        for(size_t r=1; r<15; r++)
        {
            //kick
            dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
            x[px_] -= dummy*K1*x[x_] *(R[r]+R[r-1])*lambda;
            x[pz_] += dummy*K1*x[z_] *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy *(R[r]+R[r-1])*lambda;
            //drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
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
                dummy*=R[r]*lambda;
                Cos=cos(dummy);
                Sin=sin(dummy);
                dummy = 2.0*Sin*( x[Sx_]*(Lambda[0]*Lambda[0]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[0]*Sin+Lambda[1]*Cos) +x[Ss_]*(Lambda[0]*Lambda[1]*Sin-Lambda[2]*Cos) );
                x[Ss_]+= 2.0*Sin*( x[Ss_]*(Lambda[1]*Lambda[1]-1.0)*Sin +dummy2*(Lambda[2]*Lambda[1]*Sin-Lambda[0]*Cos) +x[Sx_]*(Lambda[0]*Lambda[1]*Sin+Lambda[2]*Cos) );
                x[Sx_]+=dummy;
            }
            else{
                dummy=2.0*R[r]*lambda*(Lambda[1]*dummy2-Lambda[2]*x[Ss_]);
                x[Ss_]+=2.0*R[r]*lambda*(Lambda[2]*x[Sx_]-Lambda[0]*dummy2);
                x[Sx_]+=dummy;
            }
            //drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
        }
        //kick
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        x[px_] -= dummy*K1*x[x_] *R[14]*lambda;
        x[pz_] += dummy*K1*x[z_] *R[14]*lambda;
        x[vt_] += dummy *R[14]*lambda;
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
