// E[1][1] = E[2][2] = Young/(1-nu*nu)
// E[[1][2]= nu*E[1][1];
// E[3][3]=1/2*E/(1+nu)
// E[1][3]=E[2][3] = 0
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Math.d/TTensor.h>

ElaLinIsoMat2D::ElaLinIsoMat2D(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  t = p->eh;
  rho = p->rho;
}

ElaLinIsoMat2D::ElaLinIsoMat2D(double _rho, double _E, double _nu, double _t)
{
  rho = _rho; nu = _nu; E = _E; t = _t;
}

void
ElaLinIsoMat2D::getStress(Tensor *_stress, Tensor &_strain, double*, double temp)
{
  SymTensor<SymTensor<double,2>,2> tm;
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  getTangentMaterial(&tm, _strain, 0);
  (*stress) = tm||strain;
}

void 
ElaLinIsoMat2D::getTangentMaterial(Tensor *_tm, Tensor &, double*)
{
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2>  *>(_tm);

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;
}

void 
ElaLinIsoMat2D::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) =(*tm)||strain;
}

void 
ElaLinIsoMat2D::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                          double *staten, double *statenp, double)
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) = (*tm)||enp;
}

void
ElaLinIsoMat2D::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                          double *staten, double *statenp, double)
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = new SymTensor<SymTensor<double,2>,2>();

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) = (*tm)||enp;

  delete tm;
}


extern LinearStrain2D<9> linStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
ElaLinIsoMat2D::getGenStrainEvaluator()
{
  return &linStrain2D;
}

extern GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
StVenantKirchhoffMat2D::getGenStrainEvaluator()
{
  return &glStrain2D;
}

