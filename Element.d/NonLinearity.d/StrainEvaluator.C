#include <Math.d/matrix.h> 
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <cstdio>

LinearStrain linearStrain;
GreenLagrangeStrain greenLagrangeStrain;
LogarithmicStrain logarithmicStrain;
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
                         const Tensor &_gradU, const Tensor &_dgradUdqk)
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

  B = dgradUdqk.symPart();
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

void 
GreenLagrangeStrain::getEBandDB(Tensor &_e, Tensor &_B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_sparse & DB = static_cast<Tensor_d2s2_Sd12s34_sparse &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(_B);
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

  int size = B.getSize();
  Tensor_d1s2_full temp2(size);
  temp2 = tgradU | dgradUdqk;

  B = dgradUdqk.symPart() + temp2.symPart();

  dgradUdqk.getSymSquare(DB);
}

void 
GreenLagrangeStrain::getE(Tensor &_e,Tensor &_gradU)
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
#endif

void
LogarithmicStrain::getEBandDB(Tensor &_e, Tensor &_B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense & DB = static_cast<Tensor_d2s2_Sd12s34_dense &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(_B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Matrix3d GradU;
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      GradU(i,j) = gradU(i,j);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  for(int i=0; i<numdofs; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<3; ++k)
        dGradUdq(i)(j,k) = dgradUdqk(i,j,k);

  Eigen::Matrix3d E;
  Eigen::Matrix3d dEdqk;
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> d2Edqkdq(numdofs);

  if(GradU.isZero()) {
    //note: the eigenvalue decomposition is apparently not differentiable in this case.
    //using LinearStrain for consistency with linear elasticity at small strains
    E = 0.5*(GradU + GradU.transpose());
    for(int k=0; k<numdofs; ++k) {
      dEdqk = 0.5*(dGradUdq[k] + dGradUdq[k].transpose());
      for(int l=0; l<3; ++l)
        for(int m=l; m<3; ++m) {
          B(k,l,m) = dEdqk(l,m);
          for(int j=0; j<=k; ++j)
            DB(j,k,l,m) = 0;
        }
    }
  }
  else {
    Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec(F.transpose()*F);

    // logarithmic (Hencky) strain, Lagrangean description
    Eigen::Vector3d lnl = dec.eigenvalues().array().log();
    E = 0.5*(dec.eigenvectors() * lnl.asDiagonal() * dec.eigenvectors().adjoint());

    // Moore-Penrose pseudo inverses of C - lambda_i*I
    double tol = 1e-15;
    Eigen::Array<Eigen::Matrix3d,3,1> mpinverse;
    for(int i=0; i<3; ++i) {
      Eigen::Vector3d singularValues = dec.eigenvalues() - Eigen::Vector3d::Constant(dec.eigenvalues()[i]);
      Eigen::Vector3d invertedSingularVals;
      for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
      mpinverse[i] = dec.eigenvectors() * invertedSingularVals.asDiagonal() * dec.eigenvectors().adjoint();
    }

    // some more precomputation...
    Eigen::Array3d linv = dec.eigenvalues().array().inverse();
    Eigen::Array3d linv2 = linv.square();

    // allocate memory for intermediate derivatives
    Eigen::Array<Eigen::Array3d,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
    Eigen::Array<Eigen::Vector3d,Eigen::Dynamic,1> dlnldq(numdofs);
    Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs), d2Cdqkdq(numdofs), dydq(numdofs), d2ydqkdq(numdofs);

    for(int k=0; k<numdofs; ++k) {

      // first derivative of C=F^T*F wrt q_k
      dCdq(k) = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

      // second derivatives of C wrt q_k,q_j
      for(int j=0; j<=k; ++j) {
        Eigen::Matrix3d d2Cdqjdqk;
        if((j-k)%3 == 0) {
          Eigen::Matrix3d tmp = dGradUdq[j].transpose()*dGradUdq[k];
          d2Cdqkdq[j] = tmp + tmp.transpose();
        }
        else d2Cdqkdq[j].setZero();
      }

      for(int i=0; i<3; ++i) {
        // first derivative of lambda_i with respect to q_k
        dldq[k][i] = dec.eigenvectors().col(i).adjoint().dot(dCdq[k]*dec.eigenvectors().col(i));

        // first derivatve of y_i with respect to q_k
        // least squares solve: -(C-lambda_i*I)^{+}(dCdqk*dec.eigenvectors().col(i))
        dydq[k].col(i) = -mpinverse[i]*(dCdq[k]*dec.eigenvectors().col(i));

        for(int j=0; j<=k; ++j) {
          // second derivatives of lambda_i with respect to q_k,q_j
          d2ldqkdq(j)[i] = dec.eigenvectors().col(i).adjoint().dot(d2Cdqkdq[j]*dec.eigenvectors().col(i) + 2*dCdq[j]*dydq[k].col(i));

          // second derivatives of y_i with repsect to q_k,q_j
          d2ydqkdq(j).col(i) = -mpinverse[i]*( (dCdq[j]-dldq[j][i]*Eigen::Matrix3d::Identity())*dydq[k].col(i)
            + (dCdq[k]-dldq[k][i]*Eigen::Matrix3d::Identity())*dydq[j].col(i) + d2Cdqkdq[j]*dec.eigenvectors().col(i) )
            - dec.eigenvectors().col(i)*dydq[j].col(i).dot(dydq[k].col(i));
        }
      }
      // first derivative of E with respect to q_k
      dlnldq[k] = linv*dldq[k];
      dEdqk = 0.5*(dydq[k]*lnl.asDiagonal()*dec.eigenvectors().adjoint()
                    + dec.eigenvectors()*dlnldq[k].asDiagonal()*dec.eigenvectors().adjoint()
                      + dec.eigenvectors()*lnl.asDiagonal()*dydq[k].adjoint());
      for(int j=0; j<=k; ++j) {
        // second derivative of E with respect to q_k,j
        Eigen::Vector3d d2lnldqkdqj = linv*d2ldqkdq(j) - linv2*dldq[j]*dldq[k];
        d2Edqkdq[j] = 0.5*(d2ydqkdq(j)*lnl.asDiagonal()*dec.eigenvectors().adjoint()
                           + dydq[k]*dlnldq[j].asDiagonal()*dec.eigenvectors().adjoint()
                             + dydq[k]*lnl.asDiagonal()*dydq[j].adjoint()
                           + dydq[j]*dlnldq[k].asDiagonal()*dec.eigenvectors().adjoint()
                             + dec.eigenvectors()*d2lnldqkdqj.asDiagonal()*dec.eigenvectors().adjoint()
                               + dec.eigenvectors()*dlnldq[k].asDiagonal()*dydq[j].adjoint()
                           + dydq[j]*lnl.asDiagonal()*dydq[k].adjoint()
                             + dec.eigenvectors()*dlnldq[j].asDiagonal()*dydq[k].adjoint()
                               + dec.eigenvectors()*lnl.asDiagonal()*d2ydqkdq(j).adjoint());
      }

      for(int l=0; l<3; ++l)
        for(int m=l; m<3; ++m) {
          B(k,l,m) = dEdqk(l,m);
          for(int j=0; j<=k; ++j)
            DB(j,k,l,m) = d2Edqkdq[j](l,m);
        }
    }
  }

  for(int i=0; i<3; ++i)
    for(int j=i; j<3; ++j)
      e(i,j) = E(i,j);
#else
  std::cerr << "Error: LogarithmicStrain requires AERO-S built with Eigen3 template library." << std::endl;
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
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      GradU(i,j) = gradU(i,j);

  Eigen::Matrix3d F = GradU + Eigen::Matrix3d::Identity();
  //Eigen::JacobiSVD<Eigen::Matrix3d,Eigen::NoQRPreconditioner> dec(F, Eigen::ComputeFullV);
  //Eigen::Matrix3d E = dec.matrixV() * dec.singularValues().array().log().matrix().asDiagonal() * dec.matrixV().adjoint();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec(F.transpose()*F);
  Eigen::Matrix3d E = 0.5*(dec.eigenvectors() * dec.eigenvalues().array().log().matrix().asDiagonal() * dec.eigenvectors().adjoint());

  for(int i=0; i<3; ++i)
    for(int j=i; j<3; ++j)
      e(i,j) = E(i,j);
#else
  std::cerr << "Error: LogarithmicStrain requires AERO-S built with Eigen3 template library." << std::endl;
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
DeformationGradient::getEBandDB(Tensor &_e, Tensor &_B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_full & B = static_cast<Tensor_d1s2_full &>(_B);
  
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

