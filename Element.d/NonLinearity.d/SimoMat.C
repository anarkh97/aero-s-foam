#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/SimoMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

SimoMat::SimoMat(double _rho, double _E, double _nu)
{
  rho = _rho;
  E = _E;
  nu = _nu;
}

NLMaterial *
SimoMat::clone() const
{
  return new SimoMat(*this);
}

void
SimoMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus
  
  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

void
SimoMat::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  using std::pow;
  using std::exp;
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(_stress);
  Tensor_d0s2 &gradU = static_cast<Tensor_d0s2 &>(_gradU);
  Eigen::Matrix3d GradU; gradU.assignTo(GradU);
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::Matrix3d C = F.transpose()*F;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec(C);
  Eigen::Array<double,3,1> eps = dec.eigenvalues().array().sqrt().log();
  double J = exp(eps[0]+eps[1]+eps[2]);
  double K = E/(3*(1-2*nu)); // bulk modulus
  double dUdJ = K*log(J)/J;
  Eigen::Array<Eigen::Matrix3d,3,1> M;
  for(int i=0; i<3; ++i) {
    M[i] = 1/dec.eigenvalues()[i]*dec.eigenvectors().col(i)*dec.eigenvectors().col(i).transpose();
  }
  S = J*dUdJ*C.inverse() + (beta[0]-J*dUdJ)*M[0]+(beta[1]-J*dUdJ)*M[1]+(beta[2]-J*dUdJ)*M[2];
#else
  std::cerr << "ERROR: SimoMat::transformStress is not implemented\n";
  S.setZero();
#endif
}

void 
SimoMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "SimoMat::getTangentMaterial is not implemented\n"; exit(-1);
}

void 
SimoMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "SimoMat::getStressAndTangentMaterial is not implemented\n"; exit(-1);
}

void 
SimoMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                   double *, double *, double, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);
  Tensor_d0s4_Ss12s34_diag *tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus

  (*tm)[0][0] = (*tm)[1][1] = (*tm)[2][2] = M;
  (*tm)[0][1] = (*tm)[1][0] = (*tm)[0][2] = (*tm)[2][0] = (*tm)[1][2] = (*tm)[2][1] = lambda;

  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

void
SimoMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                   double *, double *, double, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus

  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

double
SimoMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double mu = E/(2*(1+nu));                  // shear modulus
  double trace = eps[0]+eps[1]+eps[2];

  return 0.5*lambda*trace*trace + mu*(eps[0]*eps[0]+eps[1]*eps[1]+eps[2]*eps[2]);
}

void
SimoMat::print(std::ostream &out) const
{
  out << "SimoElastic " << rho << " " << E << " " << nu;
}

void
SimoMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2);
  c[0] = E;
  c[1] = nu;
}

extern LogarithmicPrincipalStretches logarithmicPrincipalStretches;

StrainEvaluator *
SimoMat::getStrainEvaluator()
{
  return &logarithmicPrincipalStretches;
}

