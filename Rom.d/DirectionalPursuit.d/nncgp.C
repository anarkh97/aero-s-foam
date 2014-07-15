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

// Sequential implementation of Non-negative Conjugate Gradient Pursuit.
// This is a non-negative constrained variant of Conjugate Gradient Pursuit, and generates identical iterates to
// Lawson & Hanson's NNLS (non-negative least squares) algorithm.
// References:
// 1. Blumensath, Thomas, and Michael E. Davies. "Gradient pursuits." Signal Processing, IEEE Transactions on 56.6 (2008): 2370-2382.
// 2. Lawson, C. L., & Hanson, R. J. (1974). Solving least squares problems (Vol. 161). Englewood Cliffs, NJ: Prentice-hall.

Eigen::VectorXd
nncgp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      double maxsze, double reltol, bool verbose, bool scaling)
{
  using namespace Eigen;

  const long int maxvec = std::min(b.rows(), (long int)(maxsze*A.cols()));
  const long int maxit = 3*A.cols();
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  VectorXd x(A.cols()), y(A.cols()), r(A.rows()), g(A.cols()), h(A.cols()), DtGDinv(maxvec), z(maxvec), g_(maxvec), a(maxvec), S(A.cols());
  x.setZero();
  y.setZero();
  r = b;
  rnorm = bnorm;
  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec), GD(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  B.setZero(); D.setZero(); BD.setZero(); GD.setZero();
  long int k = 0;
  std::vector<long int> indices;

  if(scaling) for(int i=0; i<A.cols(); ++i) S[i] = 1/A.col(i).norm();
  else S.setOnes();

  int iter = 0; // number of iterations
  while(true) {

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k == maxvec || iter >= maxit) break; // XXX check other termination conditions

    g = S.asDiagonal()*(A.transpose()*r);
    long int i;
    h = g; for(long int j=0; j<k; ++j) h[indices[j]] = std::numeric_limits<double>::min(); // make sure the index has not already been selected?
    double gi = h.maxCoeff(&i); // note it is not maxAbs, this is different to CGP/OMP
    if(gi <= 0) break;
    B.col(k) = S[i]*A.col(i);
    indices.push_back(i);
    // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
    GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);

    for(long int j=0; j<k+1; ++j) g_[j] = g[indices[j]];
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
    d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k)*DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1));

    c = B.leftCols(k+1)*d.head(k+1);
    GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k];
    y = x; for(long int j=0; j<k+1; ++j) y[indices[j]] += a[k]*d[j];
    r -= a[k]*c;
    k++;
    while(true) {
      iter++;
      long int j;
      double yj;
      if((yj = y.minCoeff(&j)) < 0) { // XXX here we could/should take more than one if they have the same value
        std::vector<long int>::iterator pos = std::find(indices.begin(), indices.end(), j);
        std::vector<long int>::iterator fol = indices.erase(pos);
        //std::cout << "removing index " << j << std::endl;

        // note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
        // this is done here by starting from the column of D corresponding to fol (because the ones before this are already G-orthogonal)
        // and then following the same procedure that is used above to construct the original basis.
        k = std::distance(indices.begin(), fol);
        VectorXd zz = D.topLeftCorner(k,k).triangularView<Upper>()*a.head(k);
        y.setZero(); for(long int j=0; j<k; ++j) y[indices[j]] += zz[j];
        r = b - BD.leftCols(k)*a.head(k);
        for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it) {
          B.col(k) = S[*it]*A.col(*it);
          g_.head(k+1) = B.leftCols(k+1).transpose()*r;
          GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);
          Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k); 
          d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k)*DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1));
          c = B.leftCols(k+1)*d.head(k+1);
          GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          for(long int j=0; j<k+1; ++j) y[indices[j]] += a[k]*d[j];
          r -= a[k]*c;
          k++;
        }
      }
      else {
        x = y;
        break;
      }
    }

    rnorm = r.norm();
  }
  return S.asDiagonal()*x;
}

#endif
