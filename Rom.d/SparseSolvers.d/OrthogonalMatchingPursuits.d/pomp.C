#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Timers.d/GetTime.h>
#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

struct double_int {
  double val;
  int rank;
};

struct long_int {
  long index;
  int sub;
};

bool operator== (const long_int& lhs, const long_int& rhs);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pnncgp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, double &dtime);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pomp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool project, double &dtime)
{
  // each A[i] is the columnwise block of the global A matrix assigned to a subdomain on this mpi process
  // each x[i] of the return value x is the corresponding row-wise block of the global solution vector
  // n is the number of columns in the global A matrix (note this is only used to define the stopping criteria)
  using namespace Eigen;
  int myrank;
#ifdef USE_MPI
  const MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &myrank);
#else
  myrank = 0;
#endif
  struct double_int s;
  struct long_int p;

  const int nsub = A.size(); // number of subdomains assigned to this mpi process
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*n));
  const long int maxit = maxite*n;
  std::vector<long int> maxlocvec(nsub); for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  Array<VectorXd,Dynamic,1> x_(nsub), y(nsub), g(nsub), h(nsub), g_(nsub), S(nsub), t(nsub);
  for(int i=0; i<nsub; ++i) {
    x_[i] = VectorXd::Zero(maxlocvec[i]);
    y[i].resize(maxlocvec[i]);
    g_[i].resize(maxlocvec[i]);
    S[i].resize(A[i].cols());
    if(scaling) for(int j=0; j<A[i].cols(); ++j) { double s = A[i].col(j).norm(); S[i][j] = (s != 0) ? 1/s : 0; }
    else S[i].setOnes();
    t[i].resize(maxlocvec[i]);
  }
  VectorXd r(m), DtGDinv(maxvec), z(maxvec), a(maxvec);
  r = b;
  rnorm = bnorm;
  info = (n < 0) ? 2 : 1;
  a.setZero();
  Array<MatrixXd,Dynamic,1> B(nsub), D(nsub), GD(nsub);
  for(int i=0; i<nsub; ++i) {
    B[i]  = MatrixXd::Zero(m,maxlocvec[i]);
    D[i]  = MatrixXd::Zero(maxlocvec[i],maxvec);
    GD[i] = MatrixXd::Zero(maxlocvec[i],maxvec);
  }
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(m,maxvec);
  MatrixXd z_loc(maxvec,nsub), c_loc(m,nsub);
  long int k = 0; // k is the dimension of the (global) set of selected indices
  std::vector<long int> l(nsub); // l[i] is the dimension of the subset of selected indices local to a subdomain.
  std::list<std::pair<int,long_int> > gindices; // global indices
  std::list<std::pair<int,long_int> > nld_indices;
  std::vector<std::vector<long int> > indices(nsub); // local indices

  Array<long int,Dynamic,1> jk(nsub);
  ArrayXd gmax(nsub), ymin(nsub), alpha(nsub);

  long int iter   = 0; // number of iterations
  long int downIt = 0;
  while(true) {

    if(myrank == 0 && verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      g[i] = S[i].asDiagonal()*(A[i].transpose()*r);
      // make sure the index has not already been selected
      h[i] = g[i]; for(long int j=0; j<l[i]; ++j) h[i][indices[i][j]] = -std::numeric_limits<double>::max();
      // also make sure near linear dependent indices are not selected
      for(std::list<std::pair<int,long_int> >::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it)
        if(it->first == myrank && it->second.sub == i) h[i][it->second.index] = -std::numeric_limits<double>::max();
      gmax[i] = (A[i].cols() > 0) ? h[i].maxCoeff(&jk[i]) : -std::numeric_limits<double>::max();
    }
    int ik; // subdomain which has the max coeff.
    s.val  = (nsub > 0) ? gmax.maxCoeff(&ik) : -std::numeric_limits<double>::max();
    s.rank = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
#endif
    if(s.rank == myrank) {
      p.index = jk[ik];
      p.sub   = ik;
    }
#ifdef USE_MPI
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
    if(s.val <= 0) break;
    if(s.rank == myrank) {
      B[ik].col(l[ik]) = S[ik][jk[ik]]*A[ik].col(jk[ik]);
      // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
      GD[ik].row(l[ik]).head(k) = B[ik].col(l[ik]).transpose()*BD.leftCols(k);
      indices[ik].push_back(jk[ik]);
      l[ik]++;
    }
    gindices.push_back(std::pair<int,long_int>(s.rank,p));

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      for(long int j=0; j<l[i]; ++j) g_[i][j] = g[i][indices[i][j]];
      z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g_[i].head(l[i]);
    }
    z.head(k) = z_loc.topRows(k).rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, z.data(), k, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

    z.head(k) = (DtGDinv.head(k).asDiagonal()*z.head(k)).eval();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
      c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
    }
    Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k];
    if(a[k] < 0) { // check for near linear dependence
      nld_indices.push_back(std::pair<int,long_int>(s.rank,p));
      gindices.pop_back();
      if(s.rank == myrank) {
        indices[ik].pop_back();
        l[ik]--;
      }
      continue;
    } else nld_indices.clear();
    r -= a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      y[i].head(l[i]) = x_[i].head(l[i]) + a[k]*d.head(l[i]);
      GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
    }
    k++;
    iter++;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      x_[i].head(l[i]) = y[i].head(l[i]);
    }

    rnorm = r.norm();
  }
  
  if(project){
    if(myrank == 0 && verbose) {
      std::cout << "*** PROJECTING ON TO SELECTED BASIS ***" << std::endl;
    }
    std::vector<Eigen::Map<Eigen::MatrixXd> > Basis(nsub, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      new (&Basis[i]) Eigen::Map<Eigen::MatrixXd>(B[i].data(),B[i].rows(),l[i]);
    }

    x_ = pnncgp(Basis, b, rnorm, n, info, maxsze, maxite, reltol, verbose, scaling, dtime);
  }

  if(myrank == 0 && verbose) std::cout.flush();

  Array<VectorXd,Dynamic,1> x(nsub);
  for(int i=0; i<nsub; ++i) {
    x[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] = S[i][indices[i][j]]*x_[i][j];
  }
  return x;
}

#endif
