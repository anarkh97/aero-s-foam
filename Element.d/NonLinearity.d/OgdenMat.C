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
  Tensor_d0s2_Ss12_diag & strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag * stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);

  double J = strain[0]*strain[1]*strain[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  (*stress)[0] = strain[1]*strain[2]*dUdJ;
  (*stress)[1] = strain[0]*strain[2]*dUdJ;
  (*stress)[2] = strain[0]*strain[1]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(strain[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(strain[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(strain[2],alpha[i]-1);
  }
}

void 
OgdenMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12_diag & strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s4_Ss12s34_diag * tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double strain00 = strain[0]*strain[0],
         strain01 = strain[0]*strain[1],
         strain02 = strain[0]*strain[2],
         strain11 = strain[1]*strain[1],
         strain12 = strain[1]*strain[2],
         strain22 = strain[2]*strain[2];
  double J = strain01*strain[2];
  double dUdJ, d2UdJ2;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  (*tm)[0][0] = strain12*strain12*d2UdJ2;
  (*tm)[1][1] = strain02*strain02*d2UdJ2;
  (*tm)[2][2] = strain01*strain01*d2UdJ2;
  (*tm)[0][1] = (*tm)[1][0] = strain12*strain02*d2UdJ2 + strain[2]*dUdJ;
  (*tm)[0][2] = (*tm)[2][0] = strain12*strain01*d2UdJ2 + strain[1]*dUdJ;
  (*tm)[1][2] = (*tm)[2][1] = strain02*strain01*d2UdJ2 + strain[0]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(strain[0],alpha[i]-2);
    (*tm)[1][1] += mu[i]*(alpha[i]-1)*pow(strain[1],alpha[i]-2);
    (*tm)[2][2] += mu[i]*(alpha[i]-1)*pow(strain[2],alpha[i]-2);
  }
}

void 
OgdenMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12_diag & strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag * stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);
  Tensor_d0s4_Ss12s34_diag * tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double strain00 = strain[0]*strain[0],
         strain01 = strain[0]*strain[1],
         strain02 = strain[0]*strain[2],
         strain11 = strain[1]*strain[1],
         strain12 = strain[1]*strain[2],
         strain22 = strain[2]*strain[2];
  double J = strain01*strain[2];
  double dUdJ, d2UdJ2;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  (*stress)[0] = strain12*dUdJ;
  (*stress)[1] = strain02*dUdJ;
  (*stress)[2] = strain01*dUdJ;
  (*tm)[0][0] = strain12*strain12*d2UdJ2;
  (*tm)[1][1] = strain02*strain02*d2UdJ2;
  (*tm)[2][2] = strain01*strain01*d2UdJ2;
  (*tm)[0][1] = (*tm)[1][0] = strain12*strain02*d2UdJ2 + strain[2]*dUdJ;
  (*tm)[0][2] = (*tm)[2][0] = strain12*strain01*d2UdJ2 + strain[1]*dUdJ;
  (*tm)[1][2] = (*tm)[2][1] = strain02*strain01*d2UdJ2 + strain[0]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(strain[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(strain[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(strain[2],alpha[i]-1);
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(strain[0],alpha[i]-2);
    (*tm)[1][1] += mu[i]*(alpha[i]-1)*pow(strain[1],alpha[i]-2);
    (*tm)[2][2] += mu[i]*(alpha[i]-1)*pow(strain[2],alpha[i]-2);
  }
}

void 
OgdenMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                    double *, double *, double, double)
{
  Tensor_d0s2_Ss12_diag &strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);
  Tensor_d0s4_Ss12s34_diag *tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double strain00 = strain[0]*strain[0],
         strain01 = strain[0]*strain[1],
         strain02 = strain[0]*strain[2],
         strain11 = strain[1]*strain[1],
         strain12 = strain[1]*strain[2],
         strain22 = strain[2]*strain[2]; 
  double J = strain01*strain[2];
  double dUdJ, d2UdJ2; 
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    d2UdJ2 = c/2*(2+1/(J*J));          break;
    case 1: dUdJ = c*(2*J-2) - d/J;  d2UdJ2 = c*2 + d/(J*J);            break;
    case 2: dUdJ = (c*log(J) - d)/J; d2UdJ2 = (c + d - c*log(J))/(J*J); break;
  }
  (*stress)[0] = strain12*dUdJ;
  (*stress)[1] = strain02*dUdJ;
  (*stress)[2] = strain01*dUdJ;
  (*tm)[0][0] = strain12*strain12*d2UdJ2;
  (*tm)[1][1] = strain02*strain02*d2UdJ2;
  (*tm)[2][2] = strain01*strain01*d2UdJ2;
  (*tm)[0][1] = (*tm)[1][0] = strain12*strain02*d2UdJ2 + strain[2]*dUdJ;
  (*tm)[0][2] = (*tm)[2][0] = strain12*strain01*d2UdJ2 + strain[1]*dUdJ;
  (*tm)[1][2] = (*tm)[2][1] = strain02*strain01*d2UdJ2 + strain[0]*dUdJ;
  for(int i=0; i<N; ++i) { 
    (*stress)[0] += mu[i]*pow(strain[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(strain[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(strain[2],alpha[i]-1);
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(strain[0],alpha[i]-2);
    (*tm)[1][1] += mu[i]*(alpha[i]-1)*pow(strain[1],alpha[i]-2);
    (*tm)[2][2] += mu[i]*(alpha[i]-1)*pow(strain[2],alpha[i]-2);
  }
}

void
OgdenMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                    double *, double *, double, double)
{
  using std::pow;
  Tensor_d0s2_Ss12_diag &strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);

  double J = strain[0]*strain[1]*strain[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  (*stress)[0] = strain[1]*strain[2]*dUdJ;
  (*stress)[1] = strain[0]*strain[2]*dUdJ;
  (*stress)[2] = strain[0]*strain[1]*dUdJ;
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(strain[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(strain[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(strain[2],alpha[i]-1);
  }
}

double
OgdenMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  using std::log;
  using std::pow;
  Tensor_d0s2_Ss12_diag &strain = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double J = strain[0]*strain[1]*strain[2];
  double U;
  switch(vol) {
    case 0: U = c/2*((J*J-1)-log(J));         break;
    case 1: U = c*(J-1)*(J-1) - d*log(J);     break;
    case 2: U = c/2*log(J)*log(J) - d*log(J); break;
  }
  double W = 0;
  for(int i=0; i<N; ++i) W += mu[i]/alpha[i]*(pow(strain[0],alpha[i])+pow(strain[1],alpha[i])+pow(strain[2],alpha[i])-3);
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

