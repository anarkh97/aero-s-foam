#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#define UNCOUPLED_FORM

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
  // TODO: this function is used for post-processing and should return the PK2 stress
  //       currently whatever stress measure is conjugate to the principal stretches
  //       is what is being returned.
  using std::pow;
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);

  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[1] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[2] = lambda[0]*lambda[1]*dUdJ;
#ifdef UNCOUPLED_FORM
  double p = pow(J,-1/3.);
  double l[3] = { p*lambda[0], p*lambda[1], p*lambda[2] };
  double beta0=0, beta1=0, beta2=0;
  for(int i=0; i<N; ++i) {
    double w[3] = { mu[i]*pow(l[0],alpha[i]), mu[i]*pow(l[1],alpha[i]), mu[i]*pow(l[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta0 += w[0]-toto;
    beta1 += w[1]-toto;
    beta2 += w[2]-toto;
  }
  (*stress)[0] += beta0/lambda[0];
  (*stress)[1] += beta1/lambda[1];
  (*stress)[2] += beta2/lambda[2];
#else
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(lambda[2],alpha[i]-1);
  }
#endif
}

void 
OgdenMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "OgdenMat::getTangentMaterial is not implemented\n"; exit(-1);
}

void 
OgdenMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "OgdenMat::getStressAndTangentMaterial is not implemented\n"; exit(-1);
}

void 
OgdenMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                    double *, double *, double, double)
{
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);
  Tensor_d0s4_Ss12s34_diag *tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

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
  (*stress)[0] = lambda12*dUdJ;
  (*stress)[1] = lambda02*dUdJ;
  (*stress)[2] = lambda01*dUdJ;
  (*tm)[0][0] = lambda12*lambda12*d2UdJ2;
  (*tm)[1][1] = lambda02*lambda02*d2UdJ2;
  (*tm)[2][2] = lambda01*lambda01*d2UdJ2;
  (*tm)[0][1] = lambda12*lambda02*d2UdJ2 + lambda[2]*dUdJ;
  (*tm)[0][2] = lambda12*lambda01*d2UdJ2 + lambda[1]*dUdJ;
  (*tm)[1][2] = lambda02*lambda01*d2UdJ2 + lambda[0]*dUdJ;
#ifdef UNCOUPLED_FORM
  double p = pow(J,-1/3.);
  double l[3] = { p*lambda[0], p*lambda[1], p*lambda[2] };
  double beta0=0, beta1=0, beta2=0;
  double gamma00=0, gamma11=0, gamma22=0, gamma01=0, gamma02=0, gamma12=0;
  for(int i=0; i<N; ++i) {
    double w[3] = { mu[i]*pow(l[0],alpha[i]), mu[i]*pow(l[1],alpha[i]), mu[i]*pow(l[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta0 += w[0]-toto;
    beta1 += w[1]-toto;
    beta2 += w[2]-toto;
    gamma00 += alpha[i]*(w[0]+toto)/3;
    gamma11 += alpha[i]*(w[1]+toto)/3;
    gamma22 += alpha[i]*(w[2]+toto)/3;
    gamma01 += alpha[i]*(-w[0]-w[1]+toto)/3;
    gamma02 += alpha[i]*(-w[0]-w[2]+toto)/3;
    gamma12 += alpha[i]*(-w[1]-w[2]+toto)/3;
  }
  (*stress)[0] += beta0/lambda[0];
  (*stress)[1] += beta1/lambda[1];
  (*stress)[2] += beta2/lambda[2];
  (*tm)[0][0] += (gamma00-beta0)/lambda00;
  (*tm)[1][1] += (gamma11-beta1)/lambda11;
  (*tm)[2][2] += (gamma22-beta2)/lambda22;
  (*tm)[0][1] += gamma01/lambda01;
  (*tm)[0][2] += gamma02/lambda02;
  (*tm)[1][2] += gamma12/lambda12;
#else
  for(int i=0; i<N; ++i) { 
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(lambda[2],alpha[i]-1);
    (*tm)[0][0] += mu[i]*(alpha[i]-1)*pow(lambda[0],alpha[i]-2);
    (*tm)[1][1] += mu[i]*(alpha[i]-1)*pow(lambda[1],alpha[i]-2);
    (*tm)[2][2] += mu[i]*(alpha[i]-1)*pow(lambda[2],alpha[i]-2);
  }
#endif
  (*tm)[1][0] = (*tm)[0][1];
  (*tm)[2][0] = (*tm)[0][2];
  (*tm)[2][1] = (*tm)[1][2];
}

void
OgdenMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                    double *, double *, double, double)
{
  using std::pow;
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);

  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ;
  switch(vol) {
    case 0: dUdJ = c/2*(2*J-1/J);    break;
    case 1: dUdJ = c*(2*J-2) - d/J;  break;
    case 2: dUdJ = (c*log(J) - d)/J; break;
  }
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[1] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[2] = lambda[0]*lambda[1]*dUdJ;
#ifdef UNCOUPLED_FORM
  double p = pow(J,-1/3.);
  double l[3] = { p*lambda[0], p*lambda[1], p*lambda[2] };
  double beta0=0, beta1=0, beta2=0;
  for(int i=0; i<N; ++i) {
    double w[3] = { mu[i]*pow(l[0],alpha[i]), mu[i]*pow(l[1],alpha[i]), mu[i]*pow(l[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta0 += w[0]-toto;
    beta1 += w[1]-toto;
    beta2 += w[2]-toto;
  }
  (*stress)[0] += beta0/lambda[0];
  (*stress)[1] += beta1/lambda[1];
  (*stress)[2] += beta2/lambda[2];
#else
  for(int i=0; i<N; ++i) {
    (*stress)[0] += mu[i]*pow(lambda[0],alpha[i]-1);
    (*stress)[1] += mu[i]*pow(lambda[1],alpha[i]-1);
    (*stress)[2] += mu[i]*pow(lambda[2],alpha[i]-1);
  }
#endif
}

double
OgdenMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  using std::log;
  using std::pow;
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double J = lambda[0]*lambda[1]*lambda[2];
  double U;
  switch(vol) {
    case 0: U = c/2*((J*J-1)-log(J));         break;
    case 1: U = c*(J-1)*(J-1) - d*log(J);     break;
    case 2: U = c/2*log(J)*log(J) - d*log(J); break;
  }
  double W = 0;
#ifdef UNCOUPLED_FORM
  double p = pow(J,-1/3.);
  double lambdatilde[3] = { p*lambda[0], p*lambda[1], p*lambda[2] };
  for(int i=0; i<N; ++i) W += mu[i]/alpha[i]*(pow(lambdatilde[0],alpha[i])+pow(lambdatilde[1],alpha[i])+pow(lambdatilde[2],alpha[i])-3);
#else
  for(int i=0; i<N; ++i) W += mu[i]/alpha[i]*(pow(lambda[0],alpha[i])+pow(lambda[1],alpha[i])+pow(lambda[2],alpha[i])-3);
#endif
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

