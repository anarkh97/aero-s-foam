#include <Math.d/matrix.h> 
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <cstdio>

LinearStrain linearStrain;
GreenLagrangeStrain greenLagrangeStrain;
LogarithmicStrain logarithmicStrain;
PrincipalStretches principalStretches;
LogarithmicPrincipalStretches logarithmicPrincipalStretches;
DeformationGradient deformationGradient;

Tensor *
LinearStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
LinearStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LinearStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LinearStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
LinearStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_null *DB = new Tensor_d2s2_Sd12s34_null(numdofs);
  return DB;
}

void 
LinearStrain::getEBandDB(Tensor &_e, Tensor & bB, Tensor &DB,
                         const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(bB);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);

  B.setZero();
  B.addSymPart(dgradUdqk);
}

void
LinearStrain::getEandB(Tensor &_e, Tensor & __B, 
                       const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);

  B.setZero();
  B.addSymPart(dgradUdqk);
}

void 
LinearStrain::getE(Tensor &_e, Tensor &_gradU)
{
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);
}

void
LinearStrain::transformStress(Tensor &_stress, Tensor &gradU, Tensor_d0s2_Ss12 &S)
{
  // do nothing: transformation is only applied for finite-strain materials
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  S = stress;
}

Tensor *
GreenLagrangeStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
GreenLagrangeStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
GreenLagrangeStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
GreenLagrangeStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
GreenLagrangeStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_sparse *DB = new Tensor_d2s2_Sd12s34_sparse(numdofs);
  return DB;
}

Tensor *
GreenLagrangeStrain::getTempInstance(int numdofs)
{
  Tensor_d1s2_full *temp2 = new Tensor_d1s2_full(numdofs);
  return temp2;
}

void 
GreenLagrangeStrain::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *_temp2)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_sparse & DB = static_cast<Tensor_d2s2_Sd12s34_sparse &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  // de/dq = 1/2(dgradUdq^t + dgradUdq + ...)
  Tensor_d0s2 temp1;
  temp1 = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp1));
  enonsym.convertToSym(e);

  Tensor_d1s2_full * temp2 = static_cast<Tensor_d1s2_full*>(_temp2);
  *temp2 = tgradU | dgradUdqk;

  B.assignSymPart(dgradUdqk, *temp2);

  dgradUdqk.getSymSquare(DB);
}

void
GreenLagrangeStrain::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *_temp2)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  // de/dq = 1/2(dgradUdq^t + dgradUdq + ...)
  Tensor_d0s2 temp1;
  temp1 = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp1));
  enonsym.convertToSym(e);

  Tensor_d1s2_full * temp2 = static_cast<Tensor_d1s2_full*>(_temp2);
  *temp2 = tgradU | dgradUdqk;

  B.assignSymPart(dgradUdqk, *temp2);
}

void 
GreenLagrangeStrain::getE(Tensor &_e, Tensor &_gradU)
{
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 temp;
  temp = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp));
  enonsym.convertToSym(e);
}

void
GreenLagrangeStrain::transformStress(Tensor &_stress, Tensor &, Tensor_d0s2_Ss12 &S)
{
  // do nothing: stress is already PK2 in this case
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  S = stress;
}

Tensor *
LogarithmicStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
LogarithmicStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LogarithmicStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LogarithmicStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
LogarithmicStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense *DB = new Tensor_d2s2_Sd12s34_dense(numdofs);
  return DB;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicStrain::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense & DB = static_cast<Tensor_d2s2_Sd12s34_dense &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d E;
  Eigen::Matrix3d dEdqk;
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> d2Edqkdq(numdofs);

  if(GradU.isZero()) {
    //note: the eigenvalue decomposition is apparently not differentiable in this case.
    //using LinearStrain for consistency with linear elasticity at small strains
    E = 0.5*(GradU + GradU.transpose());
    for(int k=0; k<numdofs; ++k) {
      dEdqk = 0.5*(dGradUdq[k] + dGradUdq[k].transpose());
      B[k] = dEdqk;
    }
    DB.setZero();
  }
  else {
    // new implementation. this should be significantly faster...
    Eigen::Matrix<double,3,3> F = GradU + Eigen::Matrix<double,3,3>::Identity();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3> > dec((F.transpose()*F).eval());

    // logarithmic (Hencky) strain, Lagrangean description
    Eigen::Matrix<double,3,1> lnl = dec.eigenvalues().array().log();
    Eigen::Matrix<double,3,3> ylnl = dec.eigenvectors()*lnl.asDiagonal();
    E = 0.5*ylnl*dec.eigenvectors().adjoint();

    // Moore-Penrose pseudo inverses of C - lambda_i*I
    double tol = std::numeric_limits<double>::epsilon();
    Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
    for(int i=0; i<3; ++i) {
      Eigen::Matrix<double,3,1> singularValues = dec.eigenvalues() - Eigen::Matrix<double,3,1>::Constant(dec.eigenvalues()[i]);
      Eigen::Matrix<double,3,1> invertedSingularVals;
      for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
      mpinverse[i] = dec.eigenvectors() * invertedSingularVals.asDiagonal() * dec.eigenvectors().adjoint();
    }

    // some more precomputation...
    Eigen::Array<double,3,1> linv = dec.eigenvalues().array().inverse();
    Eigen::Array<double,3,1> linv2 = linv.square();

    // allocate memory for intermediate derivatives
    Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
    Eigen::Array<Eigen::Matrix<double,3,1>,Eigen::Dynamic,1> dlnldq(numdofs);
    Eigen::Array<Eigen::Matrix<double,3,3>,Eigen::Dynamic,1> dCdq(numdofs), d2Cdqkdq(numdofs), dydq(numdofs), d2ydqkdq(numdofs);
    Eigen::Matrix<double,3,3> tmp1,tmp3,toto;
    Eigen::Array<Eigen::Matrix<double,3,3>,Eigen::Dynamic,1> tmp2(numdofs);

    for(int k=0; k<numdofs; ++k) {

      // first derivative of C=F^T*F wrt q_k
      dCdq(k) = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

      // second derivatives of C wrt q_k,q_j
      for(int j=0; j<=k; ++j) {
        if((j-k)%3 == 0) {
          Eigen::Matrix<double,3,3> tmp = dGradUdq[j].transpose()*dGradUdq[k];
          d2Cdqkdq[j] = tmp + tmp.transpose();
        }
        //else d2Cdqkdq[j].setZero();
      }

      Eigen::Matrix<double,3,3> dCdqky = dCdq[k]*dec.eigenvectors();

      // first derivative of lambda with respect to q_k
      dldq[k] = (dec.eigenvectors().adjoint()*dCdqky).diagonal();

      // first derivatve of y with respect to q_k
      for(int i=0; i<3; ++i) dydq[k].col(i) = -mpinverse[i]*dCdqky.col(i);

      for(int j=0; j<=k; ++j) {
        Eigen::Matrix<double,3,3> dCdqjdydqk = dCdq[j]*dydq[k];
        if((j-k)%3 == 0) {
          Eigen::Matrix<double,3,3> d2Cdqkdqjy = d2Cdqkdq[j]*dec.eigenvectors();

          // second derivatives of lambda with respect to q_k, q_j for non-zero d2Cdqkdq[j]
          d2ldqkdq(j) = (dec.eigenvectors().adjoint()*(d2Cdqkdqjy + 2*dCdqjdydqk)).diagonal();

          // second derivatives of y with respect to q_k, q_j 
          toto = dCdqjdydqk + dCdq[k]*dydq[j] + d2Cdqkdqjy
                 - dydq[k]*dldq[j].matrix().asDiagonal() - dydq[j]*dldq[k].matrix().asDiagonal();
        }
        else {
          // second derivatives of lambda with respect to q_k, q_j for zero d2Cdqkdq[j]
          d2ldqkdq(j) = (dec.eigenvectors().adjoint()*(2*dCdqjdydqk)).diagonal();

          toto = dCdqjdydqk + dCdq[k]*dydq[j] 
                 - dydq[k]*dldq[j].matrix().asDiagonal() - dydq[j]*dldq[k].matrix().asDiagonal();
        }
        // second derivatives of y with respect to q_k, q_j
        d2ydqkdq(j) = -dec.eigenvectors()*(dydq[j].adjoint()*dydq[k]).diagonal().asDiagonal();
        for(int i=0; i<3; ++i)
          d2ydqkdq(j).col(i) -= mpinverse[i]*toto.col(i);
      }
      // first derivative of E with respect to q_k
      dlnldq[k] = linv*dldq[k];
      tmp1 = dydq[k]*lnl.asDiagonal();
      tmp2[k] = dlnldq[k].asDiagonal()*dec.eigenvectors().adjoint();
      dEdqk.triangularView<Eigen::Upper>() = 0.5*(dydq[k]*ylnl.adjoint()
                    + dec.eigenvectors()*dlnldq[k].asDiagonal()*dec.eigenvectors().adjoint()
                      + ylnl*dydq[k].adjoint());
      for(int j=0; j<=k; ++j) {
        // second derivative of E with respect to q_k,j
        Eigen::Matrix<double,3,1> d2lnldqkdqj = linv*d2ldqkdq(j) - linv2*dldq[j]*dldq[k];
        tmp3 = d2ydqkdq(j)*ylnl.adjoint() + dydq[k]*tmp2[j] + dydq[j]*(tmp1.adjoint() + tmp2[k]);
        d2Edqkdq[j].triangularView<Eigen::Upper>() = 0.5*(tmp3 + tmp3.adjoint()
                             + dec.eigenvectors()*d2lnldqkdqj.asDiagonal()*dec.eigenvectors().adjoint());
      }

      B[k] = dEdqk;
      for(int j=0; j<=k; ++j) DB[j*(2*numdofs-j-1)/2+k] = d2Edqkdq[j];
    }
  }

  e = E;

#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicStrain::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d E;
  Eigen::Matrix3d dEdqk;

  if(GradU.isZero()) {
    //note: the eigenvalue decomposition is apparently not differentiable in this case.
    //using LinearStrain for consistency with linear elasticity at small strains
    E = 0.5*(GradU + GradU.transpose());
    for(int k=0; k<numdofs; ++k) {
      dEdqk = 0.5*(dGradUdq[k] + dGradUdq[k].transpose());
      B[k] = dEdqk;
    }
  }
  else {
    // new implementation. this should be significantly faster...
    Eigen::Matrix<double,3,3> F = GradU + Eigen::Matrix<double,3,3>::Identity();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3> > dec((F.transpose()*F).eval());

    // logarithmic (Hencky) strain, Lagrangean description
    Eigen::Matrix<double,3,1> lnl = dec.eigenvalues().array().log();
    Eigen::Matrix<double,3,3> ylnl = dec.eigenvectors()*lnl.asDiagonal();
    E = 0.5*ylnl*dec.eigenvectors().adjoint();

    // Moore-Penrose pseudo inverses of C - lambda_i*I
    double tol = std::numeric_limits<double>::epsilon();
    Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
    for(int i=0; i<3; ++i) {
      Eigen::Matrix<double,3,1> singularValues = dec.eigenvalues() - Eigen::Matrix<double,3,1>::Constant(dec.eigenvalues()[i]);
      Eigen::Matrix<double,3,1> invertedSingularVals;
      for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
      mpinverse[i] = dec.eigenvectors() * invertedSingularVals.asDiagonal() * dec.eigenvectors().adjoint();
    }

    // some more precomputation...
    Eigen::Array<double,3,1> linv = dec.eigenvalues().array().inverse();

    // allocate memory for intermediate derivatives
    Eigen::Array<double,3,1> dldqk, dlnldqk;
    Eigen::Matrix3d dCdqk, dCdqky, dydqk;

    for(int k=0; k<numdofs; ++k) {

      // first derivative of C=F^T*F wrt q_k
      dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

      dCdqky = dCdqk*dec.eigenvectors();

      // first derivative of lambda with respect to q_k
      dldqk = (dec.eigenvectors().adjoint()*dCdqky).diagonal();

      // first derivatve of y with respect to q_k
      for(int i=0; i<3; ++i) dydqk.col(i) = -mpinverse[i]*dCdqky.col(i);

      // first derivative of E with respect to q_k
      dlnldqk = linv*dldqk;
      dEdqk.triangularView<Eigen::Upper>() = 0.5*(dydqk*ylnl.adjoint()
                    + dec.eigenvectors()*dlnldqk.matrix().asDiagonal()*dec.eigenvectors().adjoint()
                      + ylnl*dydqk.adjoint());

      B[k] = dEdqk;
    }
  }

  e = E;

#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicStrain::getE(Tensor &_e, Tensor &_gradU)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  // in this case e = matrix logarithm of the right stretch tensor, or alternatively 1/2 ln(F^T*F)
  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  //Eigen::JacobiSVD<Eigen::Matrix3d,Eigen::NoQRPreconditioner> dec(F, Eigen::ComputeFullV);
  //Eigen::Matrix3d E = dec.matrixV() * dec.singularValues().array().log().matrix().asDiagonal() * dec.matrixV().adjoint();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  Eigen::Matrix3d E = 0.5*(dec.eigenvectors() * dec.eigenvalues().array().log().matrix().asDiagonal() * dec.eigenvectors().adjoint());

  e = E;
#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicStrain::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  Tensor_d0s2 &gradU = static_cast<Tensor_d0s2 &>(_gradU);
  Eigen::Matrix3d GradU; gradU.assignTo(GradU);
  Eigen::Matrix3d T; stress.assignTo(T); // rotated Kirchhoff stress, T = R^{t}*tau*R --> P = R*T*R^{t}*F^{-t}
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity(); // deformation gradient
  Eigen::JacobiSVD<Eigen::Matrix3d,Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeFullV);
  Eigen::Matrix3d R = svd.matrixU()*svd.matrixV().adjoint(); // right stretch tensor
  Eigen::Matrix3d P = R*T*R.transpose()*F.transpose().inverse(); // first Piola-Kirchhoff stress tensor
  S = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
PrincipalStretches::getTMInstance()
{
  Tensor_d0s4_Ss12s34_diag *s = new Tensor_d0s4_Ss12s34_diag;
  return s;
}

Tensor *
PrincipalStretches::getStressInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
PrincipalStretches::getStrainInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
PrincipalStretches::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23_diag *B = new Tensor_d1s2_Ss23_diag(numdofs);
  return B;
}

Tensor *
PrincipalStretches::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense_diag *DB = new Tensor_d2s2_Sd12s34_dense_diag(numdofs);
  return DB;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
PrincipalStretches::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense_diag & DB = static_cast<Tensor_d2s2_Sd12s34_dense_diag &>(_DB);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);
  if(GradU.isZero()) GradU = std::numeric_limits<double>::epsilon()*Eigen::Matrix3d::Random(); // XXX

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());

  Eigen::Array<double,3,1> l = dec.eigenvalues();
  // principal stretches
  Eigen::Array<double,3,1> sqrtl = l.sqrt();
  e = sqrtl;

  // Moore-Penrose pseudo inverses of C - lambda_i*I
  double tol = std::numeric_limits<double>::epsilon();
  Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
  for(int i=0; i<3; ++i) {
    Eigen::Vector3d singularValues = dec.eigenvalues() - Eigen::Matrix<double,3,1>::Constant(dec.eigenvalues()[i]);
    Eigen::Vector3d invertedSingularVals;
    for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
    mpinverse[i] = dec.eigenvectors() * invertedSingularVals.asDiagonal() * dec.eigenvectors().adjoint();
  }

  // some more precomputation...
  Eigen::Array<double,3,1> dsqrtldl = 0.5/sqrtl;
  Eigen::Array<double,3,1> d2sqrtldl2 = -0.5*dsqrtldl/l;

  // allocate memory for intermediate derivatives
  Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs);
  Eigen::Matrix3d dCdqky, d2Cdqkdqj, dydqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C=F^T*F wrt q_k
    dCdq[k] = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of eigenvalues of C with respect to q_k
    dCdqky = dCdq[k]*dec.eigenvectors();
    dldq[k] = (dec.eigenvectors().adjoint()*dCdqky).diagonal();

    // first derivative of eigenvectors of C with respect to q_k
    for(int i=0; i<3; ++i) dydqk.col(i) = -mpinverse[i]*dCdqky.col(i);

    for(int j=0; j<=k; ++j) {
      if((j-k)%3 == 0) {
        // second derivatives of C wrt q_k,q_j
        d2Cdqkdqj = dGradUdq[j].transpose()*dGradUdq[k] + (dGradUdq[j].transpose()*dGradUdq[k]).transpose();

        // second derivatives of eigenvalues of C with respect to q_k, q_j for non-zero d2Cdqkdqj
        d2ldqkdq[j] = (dec.eigenvectors().adjoint()*(d2Cdqkdqj*dec.eigenvectors() + 2*dCdq[j]*dydqk)).diagonal();
      }
      else {
        // second derivatives of eigenvalues of C with respect to q_k, q_j for zero d2Cdqkdqj
        d2ldqkdq[j] = (dec.eigenvectors().adjoint()*(2*dCdq[j]*dydqk)).diagonal();
      }
    }

    // first derivative of principal stretches with respect to q_k
    B[k] = (dsqrtldl*dldq[k]).eval();
    for(int j=0; j<=k; ++j) {
      // second derivative of principal stretches with respect to q_k,j
      DB[j*(2*numdofs-j-1)/2+k] = (dsqrtldl*d2ldqkdq[j] + d2sqrtldl2*dldq[j]*dldq[k]).eval();
    }
  }

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
PrincipalStretches::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());

  // principal stretches
  Eigen::Array<double,3,1> sqrtl = dec.eigenvalues().array().sqrt();
  e = sqrtl;

  // some more precomputation...
  Eigen::Array<double,3,1> dsqrtldl = 0.5/sqrtl;

  // allocate memory for intermediate derivatives
  Eigen::Array<double,3,1> dldqk;
  Eigen::Matrix3d dCdqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C=F^T*F wrt q_k
    dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of eigenvalues of C with respect to q_k
    dldqk = (dec.eigenvectors().adjoint()*dCdqk*dec.eigenvectors()).diagonal();

    // first derivative of principal stretches with respect to q_k
    B[k] = (dsqrtldl*dldqk).eval();
  }

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
PrincipalStretches::getE(Tensor &_e, Tensor &_gradU)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // in this case e = principal stretches, i.e. sqrts of the eigenvalues of (F^T*F)
  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  Eigen::Array3d sqrtl = dec.eigenvalues().array().sqrt();
  e = sqrtl;

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
PrincipalStretches::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Tensor_d0s2_Ss12_diag &stress = static_cast<Tensor_d0s2_Ss12_diag &>(_stress);
  Tensor_d0s2 &gradU = static_cast<Tensor_d0s2 &>(_gradU);
  Eigen::Matrix3d GradU; gradU.assignTo(GradU);
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  Eigen::Array3d sqrtl = dec.eigenvalues().array().sqrt();
  S = stress[0]/sqrtl[0]*dec.eigenvectors().col(0)*dec.eigenvectors().col(0).transpose()
    + stress[1]/sqrtl[1]*dec.eigenvectors().col(1)*dec.eigenvectors().col(1).transpose()
    + stress[2]/sqrtl[2]*dec.eigenvectors().col(2)*dec.eigenvectors().col(2).transpose();
#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
LogarithmicPrincipalStretches::getTMInstance()
{
  Tensor_d0s4_Ss12s34_diag *s = new Tensor_d0s4_Ss12s34_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getStressInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getStrainInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23_diag *B = new Tensor_d1s2_Ss23_diag(numdofs);
  return B;
}

Tensor *
LogarithmicPrincipalStretches::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense_diag *DB = new Tensor_d2s2_Sd12s34_dense_diag(numdofs);
  return DB;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicPrincipalStretches::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense_diag & DB = static_cast<Tensor_d2s2_Sd12s34_dense_diag &>(_DB);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);
  if(GradU.isZero()) GradU = std::numeric_limits<double>::epsilon()*Eigen::Matrix3d::Random(); // XXX

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());

  Eigen::Array<double,3,1> l = dec.eigenvalues();
  // principal stretches
  Eigen::Array<double,3,1> sqrtl = l.sqrt();
  // logarithmic principal stretches
  e = sqrtl.log().eval();

  // Moore-Penrose pseudo inverses of C - lambda_i*I
  double tol = std::numeric_limits<double>::epsilon();
  Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
  for(int i=0; i<3; ++i) {
    Eigen::Vector3d singularValues = dec.eigenvalues() - Eigen::Matrix<double,3,1>::Constant(dec.eigenvalues()[i]);
    Eigen::Vector3d invertedSingularVals;
    for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
    mpinverse[i] = dec.eigenvectors() * invertedSingularVals.asDiagonal() * dec.eigenvectors().adjoint();
  }

  // some more precomputation...
  Eigen::Array<double,3,1> dedl = 0.5/l;
  Eigen::Array<double,3,1> d2edl2 = -dedl/l;

  // allocate memory for intermediate derivatives
  Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs);
  Eigen::Matrix3d dCdqky, dydqk, d2Cdqkdqj;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C=F^T*F wrt q_k
    dCdq[k] = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of eigenvalues of C with respect to q_k
    dCdqky = dCdq[k]*dec.eigenvectors();
    dldq[k] = (dec.eigenvectors().adjoint()*dCdqky).diagonal();

    // first derivative of eigenvectors of C with respect to q_k
    for(int i=0; i<3; ++i) dydqk.col(i) = -mpinverse[i]*dCdqky.col(i);

    for(int j=0; j<=k; ++j) {
      if((j-k)%3 == 0) {
        // second derivatives of C wrt q_k,q_j
        d2Cdqkdqj = dGradUdq[j].transpose()*dGradUdq[k] + (dGradUdq[j].transpose()*dGradUdq[k]).transpose();

        // second derivatives of eigenvalues of C with respect to q_k, q_j for non-zero d2Cdqkdqj
        d2ldqkdq[j] = (dec.eigenvectors().adjoint()*(d2Cdqkdqj*dec.eigenvectors() + 2*dCdq[j]*dydqk)).diagonal();
      }
      else {
        // second derivatives of eigenvalues of C with respect to q_k, q_j for zero d2Cdqkdqj
        d2ldqkdq[j] = (dec.eigenvectors().adjoint()*(2*dCdq[j]*dydqk)).diagonal();
      }
    }

    // first derivative of logarithmic principal stretches with respect to q_k
    B[k] = (dedl*dldq[k]).eval();
    for(int j=0; j<=k; ++j) {
      // second derivative of logarithmic principal stretches with respect to q_k,j
      DB[j*(2*numdofs-j-1)/2+k] = (dedl*d2ldqkdq[j] + d2edl2*dldq[j]*dldq[k]).eval();
    }
  }

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicPrincipalStretches::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  
  Eigen::Array<double,3,1> l = dec.eigenvalues();
  // principal stretches
  Eigen::Array<double,3,1> sqrtl = l.sqrt();
  // logarithmic principal stretches
  e = sqrtl.log().eval();

  // some more precomputation...
  Eigen::Array<double,3,1> dedl = 0.5/l;

  // allocate memory for intermediate derivatives
  Eigen::Array<double,3,1> dldqk;
  Eigen::Matrix3d dCdqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C=F^T*F wrt q_k
    dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of eigenvalues of C with respect to q_k
    dldqk = (dec.eigenvectors().adjoint()*dCdqk*dec.eigenvectors()).diagonal();

    // first derivative of logarithmic principal stretches with respect to q_k
    B[k] = (dedl*dldqk).eval();
  }

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicPrincipalStretches::getE(Tensor &_e, Tensor &_gradU)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // in this case e = logarithmic principal stretches, i.e. log of sqrts of the eigenvalues of (F^T*F)
  Eigen::Matrix3d GradU;
  gradU.assignTo(GradU);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  Eigen::Array3d sqrtl = dec.eigenvalues().array().sqrt();
  e = sqrtl.log().eval();

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicPrincipalStretches::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(_stress);
  Tensor_d0s2 &gradU = static_cast<Tensor_d0s2 &>(_gradU);
  Eigen::Matrix3d GradU; gradU.assignTo(GradU);
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  S = beta[0]/dec.eigenvalues()[0]*dec.eigenvectors().col(0)*dec.eigenvectors().col(0).transpose()
    + beta[1]/dec.eigenvalues()[1]*dec.eigenvectors().col(1)*dec.eigenvectors().col(1).transpose()
    + beta[2]/dec.eigenvalues()[2]*dec.eigenvectors().col(2)*dec.eigenvectors().col(2).transpose();
#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
DeformationGradient::getTMInstance()
{
  Tensor_d0s4 *s = new Tensor_d0s4;
  return s;
}

Tensor *
DeformationGradient::getStressInstance()
{
  Tensor_d0s2 *s = new Tensor_d0s2;
  return s;
}

Tensor *
DeformationGradient::getStrainInstance()
{
  Tensor_d0s2 *s = new Tensor_d0s2;
  return s;
}

Tensor *
DeformationGradient::getBInstance(int numdofs)
{
  Tensor_d1s2_full *B = new Tensor_d1s2_full(numdofs);
  return B;
}

Tensor *
DeformationGradient::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_null *DB = new Tensor_d2s2_Sd12s34_null(numdofs); // XXXX
  return DB;
}

void 
DeformationGradient::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_full & B = static_cast<Tensor_d1s2_full &>(__B);
  
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);

  // e = gradU + identity (nonsymmetric)
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  e = gradU + identity;

  B = dgradUdqk;
}

void 
DeformationGradient::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_full & B = static_cast<Tensor_d1s2_full &>(__B);
  
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);

  // e = gradU + identity (nonsymmetric)
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  e = gradU + identity;

  B = dgradUdqk;
}

void
DeformationGradient::getE(Tensor &_e,Tensor &_gradU)
{
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  // in this case e = gradU + identity (nonsymmetric)
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  e = gradU + identity;
}

void
DeformationGradient::transformStress(Tensor &_stress, Tensor &_gradU, Tensor_d0s2_Ss12 &S)
{
  Tensor_d0s2 &stress = static_cast<Tensor_d0s2 &>(_stress);
  Tensor_d0s2 &gradU  = static_cast<Tensor_d0s2 &>(_gradU);

#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > GradU(&gradU[0]), P(&stress[0]);
  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  S = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  std::cerr << "ERROR: DeformationGradient::transformStress is not implemented\n";
  S.setZero();
#endif
}

