#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template<int _N, int _M>
OgdenMat::OgdenMat(double _rho, double (&_mu)[_N], double (&_alpha)[_N], double (&_K)[_M])
{
  rho = _rho;
  N = std::min(_N,9);
  for(int i=0; i<N; ++i) { mu[i] = _mu[i]; alpha[i] = _alpha[i]; }
  M = std::min(_M,9);
  for(int i=0; i<M; ++i) K[i] = _K[i];
}

NLMaterial *
OgdenMat::clone() const
{
  return new OgdenMat(*this);
}

void
OgdenMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  using std::pow;
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag *stress = static_cast<Tensor_d0s2_Ss12_diag *>(_stress);

  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ = 0;
  for(int i=0; i<M; ++i) dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[1] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[2] = lambda[0]*lambda[1]*dUdJ;

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
}

void
OgdenMat::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  using std::pow;
  Tensor_d0s2_Ss12_diag &stress = static_cast<Tensor_d0s2_Ss12_diag &>(_stress);
  Tensor_d0s2 &gradU = static_cast<Tensor_d0s2 &>(_gradU);
  Eigen::Matrix3d GradU; gradU.assignTo(GradU);
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::Matrix3d C = F.transpose()*F;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec(C);
  Eigen::Array<double,3,1> lambda = dec.eigenvalues().array().sqrt();
  double J = lambda[0]*lambda[1]*lambda[2];
  double dUdJ = 0;
  for(int i=0; i<M; ++i) dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
  Eigen::Array<Eigen::Matrix3d,3,1> M;
  Eigen::Array<double,3,1> beta;
  for(int i=0; i<3; ++i) {
    M[i] = 1/dec.eigenvalues()[i]*dec.eigenvectors().col(i)*dec.eigenvectors().col(i).transpose();
    beta[i] = (stress[i]*lambda[i]-J*dUdJ);
  }
  S = J*dUdJ*C.inverse() + beta[0]*M[0]+beta[1]*M[1]+beta[2]*M[2];
#else
  std::cerr << "ERROR: OgdenMat::transformStress is not implemented\n";
  S.setZero();
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
  double dUdJ = 0, d2UdJ2 = 0; 
  for(int i=0; i<M; ++i) {
    dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
    d2UdJ2 += K[i]*(i+1)*(2*i+1)*pow(J-1,2*i);
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
  double dUdJ = 0;
  for(int i=0; i<M; ++i) dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
  (*stress)[0] = lambda[1]*lambda[2]*dUdJ;
  (*stress)[1] = lambda[0]*lambda[2]*dUdJ;
  (*stress)[2] = lambda[0]*lambda[1]*dUdJ;

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
}

double
OgdenMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  using std::log;
  using std::pow;
  Tensor_d0s2_Ss12_diag &lambda = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double J = lambda[0]*lambda[1]*lambda[2];
  double U = 0;
  for(int i=0; i<M; ++i) U += K[i]/2*pow(J-1,2*i+2);

  double W = 0;
  double p = pow(J,-1/3.);
  double lambdatilde[3] = { p*lambda[0], p*lambda[1], p*lambda[2] };
  for(int i=0; i<N; ++i) W += mu[i]/alpha[i]*(pow(lambdatilde[0],alpha[i])+pow(lambdatilde[1],alpha[i])+pow(lambdatilde[2],alpha[i])-3);

  return W+U;
}

void
OgdenMat::print(std::ostream &out) const
{
  out << "Ogden " << rho;
  for(int i=0; i<N; ++i) out << " " << mu[i];
  for(int i=0; i<N; ++i) out << " " << alpha[i];
  for(int i=0; i<M; ++i) out << " " << K[i];
}

void
OgdenMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2*N);
  for(int i=0; i<N; ++i) { c[i] = mu[i]; c[N+i] = alpha[i]; }
}

extern PrincipalStretches principalStretches;

StrainEvaluator *
OgdenMat::getStrainEvaluator()
{
  return &principalStretches;
}

template OgdenMat::OgdenMat(double, double (&)[1], double (&)[1], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[2], double (&)[2], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[3], double (&)[3], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[4], double (&)[4], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[5], double (&)[5], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[6], double (&)[6], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[7], double (&)[7], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[8], double (&)[8], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[9], double (&)[9], double (&)[1]);

template OgdenMat::OgdenMat(double, double (&)[2], double (&)[2], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[3], double (&)[3], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[4], double (&)[4], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[5], double (&)[5], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[6], double (&)[6], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[7], double (&)[7], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[8], double (&)[8], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[9], double (&)[9], double (&)[2]);
