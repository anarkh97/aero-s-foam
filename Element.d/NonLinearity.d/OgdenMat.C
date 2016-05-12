#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

OgdenMat::OgdenMat(double _rho, double _mu1, double _alpha1,
                   double _c, double _d, int _vol)
{
  rho = _rho;
  N = 1;
  mu[0] = _mu1;
  mu[1] = 0;
  mu[2] = 0;
  alpha[0] = _alpha1;
  alpha[1] = 0;
  alpha[2] = 0;
  c = _c;
  d = _d;
  vol = _vol;
}

OgdenMat::OgdenMat(double _rho, double _mu1, double _mu2, double _alpha1, double _alpha2,
                   double _c, double _d, int _vol)
{
  rho = _rho;
  N = 2;
  mu[0] = _mu1;
  mu[1] = _mu2;
  mu[2] = 0;
  alpha[0] = _alpha1;
  alpha[1] = _alpha2;
  alpha[2] = 0;
  c = _c;
  d = _d;
  vol = _vol;
}

OgdenMat::OgdenMat(double _rho, double _mu1, double _mu2, double _mu3, double _alpha1, double _alpha2, double _alpha3,
                   double _c, double _d, int _vol)
{
  rho = _rho;
  N = 3;
  mu[0] = _mu1;
  mu[1] = _mu2;
  mu[2] = _mu3;
  alpha[0] = _alpha1;
  alpha[1] = _alpha2;
  alpha[2] = _alpha3;
  c = _c;
  d = _d;
  vol = _vol;
}

NLMaterial *
OgdenMat::clone() const
{
  return new OgdenMat(*this);
}

void
OgdenMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  // Note: this function is used for post-processing and should return the PK2 stress
  //       currently whatever stress measure is conjugate to the principal stretches
  //       is what is being returned.
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  double lambda[3] = {strain[0], strain[3], strain[5]};
  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  stress->setZero();
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[3] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[5] = lambda[0]*lambda[1]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[3] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[5] += mu[i]*pow(lambda[2],alpha[i]-1);
  }
}

void 
OgdenMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda[3] = {strain[0], strain[3], strain[5]};
  double lambda00 = lambda[0]*lambda[0],
         lambda01 = lambda[0]*lambda[1],
         lambda02 = lambda[0]*lambda[2],
         lambda11 = lambda[1]*lambda[1],
         lambda12 = lambda[1]*lambda[2],
         lambda22 = lambda[2]*lambda[2];
  double J = lambda01*lambda[2];
  double dUdJ, d2UdJ2;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  tm->setZero();
  (*tm)[0][0] = lambda12*lambda12*d2UdJ2;
  (*tm)[3][3] = lambda02*lambda02*d2UdJ2;
  (*tm)[5][5] = lambda01*lambda01*d2UdJ2;
  (*tm)[0][3] = (*tm)[3][0] = lambda12*lambda02*d2UdJ2 + lambda[2]*dUdJ;
  (*tm)[0][5] = (*tm)[5][0] = lambda12*lambda01*d2UdJ2 + lambda[1]*dUdJ;
  (*tm)[3][5] = (*tm)[5][3] = lambda02*lambda01*d2UdJ2 + lambda[0]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(lambda[0],alpha[i]-2);
    (*tm)[3][3] += mu[i]*(alpha[i]-1)*pow(lambda[1],alpha[i]-2);
    (*tm)[5][5] += mu[i]*(alpha[i]-1)*pow(lambda[2],alpha[i]-2);
  }
}

void 
OgdenMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda[3] = {strain[0], strain[3], strain[5]};
  double lambda00 = lambda[0]*lambda[0],
         lambda01 = lambda[0]*lambda[1],
         lambda02 = lambda[0]*lambda[2],
         lambda11 = lambda[1]*lambda[1],
         lambda12 = lambda[1]*lambda[2],
         lambda22 = lambda[2]*lambda[2];
  double J = lambda01*lambda[2];
  double dUdJ, d2UdJ2;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  stress->setZero();
  tm->setZero();
  (*stress)[0] = lambda12*dUdJ;
  (*stress)[3] = lambda02*dUdJ;
  (*stress)[5] = lambda01*dUdJ;
  (*tm)[0][0] = lambda12*lambda12*d2UdJ2;
  (*tm)[3][3] = lambda02*lambda02*d2UdJ2;
  (*tm)[5][5] = lambda01*lambda01*d2UdJ2;
  (*tm)[0][3] = (*tm)[3][0] = lambda12*lambda02*d2UdJ2 + lambda[2]*dUdJ;
  (*tm)[0][5] = (*tm)[5][0] = lambda12*lambda01*d2UdJ2 + lambda[1]*dUdJ;
  (*tm)[3][5] = (*tm)[5][3] = lambda02*lambda01*d2UdJ2 + lambda[0]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[3] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[5] += mu[i]*pow(lambda[2],alpha[i]-1);
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(lambda[0],alpha[i]-2);
    (*tm)[3][3] += mu[i]*(alpha[i]-1)*pow(lambda[1],alpha[i]-2);
    (*tm)[5][5] += mu[i]*(alpha[i]-1)*pow(lambda[2],alpha[i]-2);
  }
}

void 
OgdenMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                    double *, double *, double, double)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda[3] = {enp[0], enp[3], enp[5]};
  double lambda00 = lambda[0]*lambda[0],
         lambda01 = lambda[0]*lambda[1],
         lambda02 = lambda[0]*lambda[2],
         lambda11 = lambda[1]*lambda[1],
         lambda12 = lambda[1]*lambda[2],
         lambda22 = lambda[2]*lambda[2]; 
  double J = lambda01*lambda[2];
  double dUdJ, d2UdJ2; 
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  stress->setZero();
  tm->setZero();
  (*stress)[0] = lambda12*dUdJ;
  (*stress)[3] = lambda02*dUdJ;
  (*stress)[5] = lambda01*dUdJ;
  (*tm)[0][0] = lambda12*lambda12*d2UdJ2;
  (*tm)[3][3] = lambda02*lambda02*d2UdJ2;
  (*tm)[5][5] = lambda01*lambda01*d2UdJ2;
  (*tm)[0][3] = (*tm)[3][0] = lambda12*lambda02*d2UdJ2 + lambda[2]*dUdJ;
  (*tm)[0][5] = (*tm)[5][0] = lambda12*lambda01*d2UdJ2 + lambda[1]*dUdJ;
  (*tm)[3][5] = (*tm)[5][3] = lambda02*lambda01*d2UdJ2 + lambda[0]*dUdJ;
  for(int i=0; i<N; ++i) { 
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[3] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[5] += mu[i]*pow(lambda[2],alpha[i]-1);
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(lambda[0],alpha[i]-2);
    (*tm)[3][3] += mu[i]*(alpha[i]-1)*pow(lambda[1],alpha[i]-2);
    (*tm)[5][5] += mu[i]*(alpha[i]-1)*pow(lambda[2],alpha[i]-2);
  }
}

void
OgdenMat::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                    double *, double *, double, double)
{
  using std::pow;
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  double lambda[3] = {enp[0], enp[3], enp[5]};
  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  stress->setZero();
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[3] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[5] = lambda[0]*lambda[1]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[3] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[5] += mu[i]*pow(lambda[2],alpha[i]-1);
  }
}

double
OgdenMat::getStrainEnergyDensity(Tensor &_enp, double *, double)
{
  using std::log;
  using std::pow;
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);

  double lambda[3] = {enp[0], enp[3], enp[5]};
  double J = lambda[0]*lambda[1]*lambda[2];
  double U;
  switch(vol) {
    case 0: U = c/2*((J*J-1)-log(J));         break;
    case 1: U = c*(J-1)*(J-1) - d*log(J);     break;
    case 2: U = c/2*log(J)*log(J) - d*log(J); break;
  }
  double W = 0;
  for(int i=0; i<N; ++i) W += mu[i]/alpha[i]*(pow(lambda[0],alpha[i])+pow(lambda[1],alpha[i])+pow(lambda[2],alpha[i])-3);
  return W+U;
}

void
OgdenMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(6);
  c[0] = mu[0];
  c[1] = mu[1];
  c[2] = mu[2];
  c[3] = alpha[0];
  c[4] = alpha[1];
  c[5] = alpha[2];
}

extern PrincipalStretches principalStretches;

StrainEvaluator *
OgdenMat::getStrainEvaluator()
{
  return &principalStretches;
}

