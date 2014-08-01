#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>

// AUTHOR: Todd Chapman, Stanford University 2014
// This function solves the Compressed Sensing problem
//
// min ||x||_0 s.t. A*x = b
//
// which is equivalent to 
//
// min ||s||_1 s.t. [A,-A]*s = b & s >= 0 
//
// and its dual problem
//
// min b^T*c s.t. [A,-A]^T*c <= 1
//
// where s_i = x_i > 0 for i = 1:N
// and   s_i = x_i < 0 for i = N+1:2*N
//
// The solution framework is known as Basis Pursuit and the algorithm is called Gradient Polytope 
// Faces Pursuit. This is a particular implementation of PFP which uses the Conjugate Gradient method to 
// solve a system of equations at each iteration whose dimensionality changes from one iteration to the 
// next, thus only matrix vector products are required and the scheme is easily parallelizable.
// This algorithm produces identical residuals to that of the standard Polytope Faces Pursuit
//
// References:
// 1. Plumbley, M & Gretsistas, I. "Gradient Polytope Faces Pursuit for Large Sparse Recovery Problems" ICASSP
// 2. Blumensath, T & Davies, M. "Gradient Pursuits" IEEE
//
// ARGUMENTS:
// A       = where A*x= b
// b       = target vector 
// rnorm   = residual norm 
// maxsze  = maximum allowable sparsity
// reltol  = stopping criteria
// scaling = flag to turn on/off unit normalization of columns of A
 

Eigen::VectorXd
gpfp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      double maxsze, double reltol, bool verbose, bool scaling)
{
  using namespace Eigen;

  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = 3*A.cols();
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  VectorXd y(maxvec), r(A.rows()), vertex(A.rows()), g1(A.cols()), g2(A.cols()), h(A.cols()), DtGDinv(maxvec), g_(maxvec), lambda(maxvec), a(maxvec), S(A.cols()), t(maxvec);
  r = b;
  vertex.setZero();
  rnorm = bnorm;
  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec), GD(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  B.setZero(); D.setZero(); BD.setZero(); GD.setZero();
  long int k = 0;
  std::vector<long int> indices;
  std::vector<int>      setKey;

  if(scaling) for(int i=0; i<A.cols(); ++i) S[i] = 1/A.col(i).norm();
  else S.setOnes();

  int iter   = 0; // number of iterations
  int downIt = 0; // number of downdates
  while(true) {

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate  = " << std::setw(9) << downIt << "   "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k == maxvec || iter >= maxit) break;

    g1.setZero();
    g2.setZero();
    int Set = 1;
    long int position  = 0;
    long int position1 = 0;
    long int position2 = 0;

    // compute step length along current face
    g1 = S.asDiagonal()*(A.transpose()*r);
    g2 = S.asDiagonal()*(A.transpose()*vertex);

    //zero out elements that are already selected
    h = g1; 
    for(long int j=0; j<k; ++j) h(indices[j]) = 0.0;

    for(long int col = 0; col != A.cols(); col++){
      double num = h(col);
      double den = g2(col);
      if(num >= 0.) {
        h(col) = num/(1.0-den);
        g2(col) = std::numeric_limits<double>::min();
      } else {
        h(col) = std::numeric_limits<double>::min();
        g2(col) = (-1.0)*num/(1.0+den);
      }
    }

    double lam  = 0.0;
    double lam1 = h.maxCoeff(&position1);
    double lam2 = g2.maxCoeff(&position2);

    if (lam1 > lam2){
      lam = lam1;
      position = position1;
    } else {
      lam = lam2;
      position = position2;
      Set = -1;
    }

    // remove maximum element from the correct active set and place
    // in correct dual set, also store appropriate column of [A,-A]
    setKey.push_back(Set);
    B.col(k) = double(Set)*S[position]*A.col(position);
    indices.push_back(position);

    // add step length for vertex estimate update to array
    lambda[k] = 1./lam;

    // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
    GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);

    for(long int j=0; j<k+1; ++j) g_[j] = g1[indices[j]];
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
    d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1))); // direction

    c = B.leftCols(k+1)*d.head(k+1);
    GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k]; // step length
    y.head(k+1) += a[k]*d.head(k+1); // candidate solution
    r -= a[k]*c; // residual
    vertex += lambda[k]*a[k]*c;
    k++;
    while(true) {
      iter++;
      int i = 0;
      double minCoeff = y.head(k).minCoeff(&i); 
      if(minCoeff < 0) {
        downIt++;
        //std::cout << " removing index " << indices[i] << std::endl;
        std::vector<long int>::iterator pos = indices.begin()+i;
        std::vector<long int>::iterator fol = indices.erase(pos);
        setKey.erase(setKey.begin()+i);

        // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
        // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
        // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
        y.head(k).setZero();
        long int k0 = k = std::distance(indices.begin(), fol);
        y.head(k) = D.topLeftCorner(k,k).triangularView<Upper>()*a.head(k);
        r = b - BD.leftCols(k)*a.head(k);
        vertex = BD.leftCols(k)*lambda.head(k).asDiagonal()*a.head(k);
        g_.head(k) = B.leftCols(k).transpose()*r;
        for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it) {
          B.col(k) = double(setKey[k])*S[*it]*A.col(*it);
          g_[k] = B.col(k).transpose()*r;
          lambda[k] = (1.0-B.col(k).dot(vertex))/(g_[k]);
          GD.row(k).head(k0) = GD.row(k+1).head(k0);
          GD.row(k).segment(k0,k-k0) = B.col(k).transpose()*BD.block(0,k0,m,k-k0);
          Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k); 
          d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1)));
          c = B.leftCols(k+1)*d.head(k+1);
          GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          y.head(k+1) += a[k]*d.head(k+1);
          r -= a[k]*c;
          vertex += lambda[k]*a[k]*c;
          k++;
          g_.head(k) -= a[k-1]*GD.col(k-1).head(k);
        }
      }
      else {
        break;
      }
    }

    rnorm = r.norm();
  }

  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  long int element  = 0;
  for(std::vector<int>::iterator it = setKey.begin(); it != setKey.end(); it++) {
    x[indices[element]] +=  S[indices[element]]*double(*it)*y[element];
    element++;
  }
  return x;
}
#endif
