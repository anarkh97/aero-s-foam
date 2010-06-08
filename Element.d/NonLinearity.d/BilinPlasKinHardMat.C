#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Math.d/mathUtility.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <limits>

BilinPlasKinHardMat::BilinPlasKinHardMat(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  Ep = p->Ep;         // Tangent modulus from the strain-stress curve
  sigE = p->sigE;     // Yield equivalent stress
}

void
BilinPlasKinHardMat::getStress(Tensor *stress, Tensor &strain, double *state)
{

}

void 
BilinPlasKinHardMat::getElasticity(Tensor *_tm)
{
  Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = lambda/nu;

  (tm)[0][0] = lambdadivnu*(1-nu);
  (tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (tm)[3][3] = lambdadivnu*(1-nu);
  (tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (tm)[5][5] = lambdadivnu*(1-nu);
  (tm)[0][3] = lambdadivnu*nu;
  (tm)[3][0] = lambdadivnu*nu;
  (tm)[0][5] = lambdadivnu*nu;
  (tm)[5][0] = lambdadivnu*nu;
  (tm)[3][5] = lambdadivnu*nu;
  (tm)[5][3] = lambdadivnu*nu;
}

void 
BilinPlasKinHardMat::getTangentMaterial(Tensor *tm, Tensor &strain, double *state)
{

}

void 
BilinPlasKinHardMat::getStressAndTangentMaterial(Tensor *stess, Tensor *tm, Tensor &strain, double *state)
{

}

void 
BilinPlasKinHardMat::updateStates(Tensor en, Tensor enp, double *state)
{

}

void
BilinPlasKinHardMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &_en, Tensor  &_enp,
                               double *staten, double *statenp, double)
{
  //////////////////////////////////////////////////////////////////////////////
  /// Simo and Hughes - Computational Inelasticity - Springer -1998- (p:124) ///
  //////////////////////////////////////////////////////////////////////////////
  // PJSA note: the presentation of this algorithm in Simo and Hughes is for
  // infintesimal strain, i.e. nonlinear material but still geometrically linear
  // I'm not sure if it is valid to simply replace the infintesimal strain with
  // the Green-Lagrange strain and use the same constitutive law.
  

  // theta == 0 corresponds to Kinematic hardening and theta == 1 to isotropic 
  // hardening
  double theta = 0.0;                          

  if(statenp == 0) {

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = lambda/nu;

    Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm); 
    Tensor_d0s2_Ss12 & enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
    Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress); 

    (*tm)[0][0] = lambdadivnu*(1-nu);
    (*tm)[1][1] = lambdadivnu*(1-2*nu)/2;
    (*tm)[2][2] = lambdadivnu*(1-2*nu)/2;
    (*tm)[3][3] = lambdadivnu*(1-nu);
    (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
    (*tm)[5][5] = lambdadivnu*(1-nu);
    (*tm)[0][3] = lambdadivnu*nu;
    (*tm)[3][0] = lambdadivnu*nu;
    (*tm)[0][5] = lambdadivnu*nu;
    (*tm)[5][0] = lambdadivnu*nu;
    (*tm)[3][5] = lambdadivnu*nu;
    (*tm)[5][3] = lambdadivnu*nu;

    (*stress) = (*tm)||enp;
  }
  else {
    //state: from 0 to 5, plastic strain ; from 6 to 11, center of the yield surface; 12 equivalent plastic strain

    Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm); 
    Tensor_d0s2_Ss12 & enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
    Tensor_d0s2_Ss12 & stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress); 

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double mu = E/(2*(1+nu));
    double bulk = lambda + (2./3)*mu;
    double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : numeric_limits<double>::max();

    Tensor_d0s2_Ss12 betan;
    Tensor_d0s2_Ss12 edevnp;
    Tensor_d0s2_Ss12 eplastn;
    Tensor_d0s2_Ss12 temp;
    Tensor_d0s2_Ss12 eplastdevn;
    Tensor_d0s2_Ss12 strialnp;
    Tensor_d0s2_Ss12 xitrialnp;

    eplastn.buildTensorOf(staten);
    betan.buildTensorOf(staten+6);

    enp.getDeviation(edevnp);
    eplastn.getDeviation(eplastdevn);

    strialnp = 2*mu*(edevnp - eplastdevn);
    xitrialnp = strialnp - betan;

    double xitrialnpnorm = sqrt(2*(xitrialnp.secondInvariant()));

    double trialyieldnp = xitrialnpnorm - sqrt(2./3)*(sigE+theta*Hprime*staten[12]);//Yield Criterion

    if (trialyieldnp <= 0) {
 
      //fprintf(stderr, "je suis dans la yield surface\n");

      for (int i=0; i<13; ++i) {
        //fprintf(stderr, "staten[%d]=%e\n", i, staten[i]);
        statenp[i] = staten[i] ;
      }

      getElasticity(_tm);

      temp = (enp - eplastn);
      (stress) = (tm) || temp;

    }
    else {
      //fprintf(stderr, "je suis en dehors de la yield surface\n");
      int i,j,k,l;

      double plastmult = trialyieldnp/(2*mu+(2./3)*Hprime); 
    
      Tensor_d0s2_Ss12 normalnp = (1./xitrialnpnorm)*xitrialnp;
      double thetanp = 1 - 2*mu*plastmult/xitrialnpnorm;
      double thetaprimenp = 1/(1+Hprime/(3*mu)) - (1 - thetanp);

      statenp[12] = staten[12] + sqrt(2./3)*plastmult;

      for (i=0; i<6; ++i) {
        statenp[i] = staten[i] + plastmult*normalnp[i];
        statenp[6+i] = staten[6+i] + sqrt(2./3)*Hprime*(1-theta)*(statenp[12]-staten[12])*normalnp[i];   
      }

      (stress) = ((bulk*3*(enp - edevnp) + strialnp) - (2*mu*plastmult*normalnp));

      for (i=0; i<3; ++i)
        for (j=i; j<3; ++j)
          for (k=0; k<3; ++k)
            for (l=k; l<3; ++l)
              tm[i*(5-i)/2+j][k*(5-k)/2+l] = (bulk*delta(i,j)*delta(k,l))
               +(2*mu*thetanp*((delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))/2.-delta(i,j)*delta(k,l)/3.))
               - 2*mu*thetaprimenp*normalnp[i*(5-i)/2+j]*normalnp[k*(5-k)/2+l];
    }
  }
}

void 
BilinPlasKinHardMat::initStates(double *st)
{
  for (int i=0; i<13; ++i)
    st[i] = 0;
}
