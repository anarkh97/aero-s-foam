#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Material.d/Material.h>

ElaLinIsoMat::ElaLinIsoMat(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  rho = p->rho; // PJSA
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double _E, double _nu)
{
  rho = _rho; nu = _nu; E = _E;
}

void
ElaLinIsoMat::getStress(Tensor *_stress,Tensor &_strain, double*)
{
  Tensor_d0s4_Ss12s34 tm;
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  getTangentMaterial(&tm, _strain, 0);
  (*stress) = tm||strain;
}

void 
ElaLinIsoMat::getTangentMaterial(Tensor *_tm,Tensor &,double*)
{
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = lambda/nu;

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
}

void 
ElaLinIsoMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = lambda/nu;

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

  (*stress) =(*tm)||strain;
}

void 
ElaLinIsoMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                        double *staten, double *statenp, double)
{
  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = lambda/nu;

  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

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
