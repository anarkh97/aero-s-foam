#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>

ElaLinIsoMat::ElaLinIsoMat(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  rho = p->rho;
  m_tm = 0;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double _E, double _nu)
{
  rho = _rho; nu = _nu; E = _E;
  m_tm = 0;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double C[6][6])
{
  rho = _rho;
  setTangentMaterial(C);
}

ElaLinIsoMat::~ElaLinIsoMat()
{
  if(m_tm) delete m_tm;
}

NLMaterial *
ElaLinIsoMat::clone() const
{
  return new ElaLinIsoMat(*this);
}

void
ElaLinIsoMat::setTangentMaterial(double C[6][6])
{
  m_tm = new Tensor_d0s4_Ss12s34();
  int index_map[6] = { 0,3,5,1,4,2 };
  for(int i=0; i<6; ++i)
    for(int j=0; j<6; ++j)
      (*m_tm)[i][j] = C[index_map[i]][index_map[j]];
}

void
ElaLinIsoMat::getStress(Tensor *_stress, Tensor &_strain, double*)
{
  Tensor_d0s4_Ss12s34 tm;
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  getTangentMaterial(&tm, _strain, 0);
  (*stress) = tm||strain;
}

void 
ElaLinIsoMat::getTangentMaterial(Tensor *_tm, Tensor &, double*)
{
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  if(m_tm) *tm = *m_tm;
  else {
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

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
}

void 
ElaLinIsoMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  if(m_tm) *tm = *m_tm;
  else {

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

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

  (*stress) =(*tm)||strain;
}

void 
ElaLinIsoMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                        double *, double *, double)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  if(m_tm) *tm = *m_tm;
  else {
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

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

  (*stress) = (*tm)||enp;
}

void
ElaLinIsoMat::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                        double *, double *, double)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  if(m_tm) {
    (*stress) = (*m_tm)||enp;
  }
  else {
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    Tensor_d0s4_Ss12s34 *tm = new Tensor_d0s4_Ss12s34();

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

    delete tm;
  }
}

double
ElaLinIsoMat::getStrainEnergyDensity(Tensor &_enp, double *)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);

  if(m_tm) {
    Tensor_d0s2_Ss12 stress;
    stress = (*m_tm)||enp;
    return 0.5*(enp[0]*stress[0] + enp[3]*stress[3] + enp[5]*stress[5]
           + 2.0*(enp[1]*stress[1] + enp[2]*stress[2] + enp[4]*stress[4]));
  }
  else {
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));

    double I1 = enp.getTrace();
    return lambda/2*I1*I1 + mu*enp.innerProduct();
  }
}

extern LinearStrain linearStrain;

StrainEvaluator *
ElaLinIsoMat::getStrainEvaluator()
{
  return &linearStrain;
}

extern GreenLagrangeStrain greenLagrangeStrain;

StrainEvaluator *
StVenantKirchhoffMat::getStrainEvaluator()
{
  return &greenLagrangeStrain;
}

NLMaterial * 
StVenantKirchhoffMat::clone() const
{
  return new StVenantKirchhoffMat(*this);
}

extern LogarithmicStrain logarithmicStrain;

StrainEvaluator *
HenckyMat::getStrainEvaluator()
{
  return &logarithmicStrain;
}

NLMaterial *
HenckyMat::clone() const
{
  return new HenckyMat(*this);
}

