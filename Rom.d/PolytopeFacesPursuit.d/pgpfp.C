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
// where s_i   = x_i > 0 for i = 1:N
// and   s_2*i = x_i < 0 for i = 1:N
//
// The solution framework is known as Basis Pursuit and the algorithm is called Gradient Polytope 
// Faces Pursuit. This is a paraticular implementation of PFP which uses the Conjugate Gradient method to 
// solve a system of equations at each iteration who's dimensionality changes from one iteration to the 
// next, thus only matrix vector products are required and the scheme is easily parallelizable.
// This algorithm produces identical residuals to that of the standard Polytope Faces Pursuit
//
// References:
// 1. Plumbley, M & Gretsistas, I. "Gradient Polytope Faces Pursuit for Large Sparse Recovery Problems" ICASSP
// 2. Blumensath, T & Davies, M. "Gradient Pursuits" IEEE
//
// AURGUMENTS:
// A       = where A*x= b and A[i] is columnwise block of global A matrix assigned to subdomain on this mpi process
// b       = target vector 
// rnorm   = residual norm 
// n       = number of columns in global A matrix (for stopping criteria only)
// maxsze  = maximum allowable sparsity
// reltol  = stopping criteria
// scaling = flag to turn on/off unit normalization of columns of A
//
// each x[i] of the return value x is the corresponding row-wise block of the global solution vector
Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pgpfp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       double maxsze, double reltol, bool verbose, bool scaling)
{
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

  const int nsub   = A.size();                               // number of subdomains assigned to this mpi process
  const long int m = b.rows();                               // number of measurements 
  const long int maxvec = std::min(m, (long int)(maxsze*n)); // maximum allowable non-zero elements
  const long int maxit  = 3*n;                               // maximum allowable iterations to ensure exit
  std::vector<long int> maxlocvec(nsub); 	             // vector of subdomain sizes 
  double bnorm  = b.norm();                                  // convergence criteria 
  double abstol = reltol*bnorm;                              // convergence criteria

  // local vectors, contain data specific to this mpi process
  Array<VectorXd,Dynamic,1> y(nsub);  // solution container
  Array<VectorXd,Dynamic,1> g1(nsub); // working space
  Array<VectorXd,Dynamic,1> g2(nsub); // working space
  Array<VectorXd,Dynamic,1> g3(nsub); // working space
  Array<VectorXd,Dynamic,1> g4(nsub); // working space
  Array<VectorXd,Dynamic,1> S(nsub);  // columns normalization container

  // global vectors, all mpi process have identical copies 
  VectorXd r(m);                      // residual
  VectorXd c(m);                      // vertex estimate
  VectorXd DtGDinv(maxvec);           // inverse of D^T*A^T*A*D (see ref. 2)
  VectorXd alpha(maxvec);             // step length container for solution and residual updates
  VectorXd lambda(maxvec);            // step length container for vertex estimate updates
  VectorXd z(maxvec);                 //

  // local matrix containers, each mpi process has a unique portion of these matrices
  Array<MatrixXd,Dynamic,1> B(nsub);  // container for selected columns of [A,-A]
  Array<MatrixXd,Dynamic,1> D(nsub);  // container for conjugate direction vectors
  Array<MatrixXd,Dynamic,1> GD(nsub); // container for D^T*G (see ref. 2)
   
  MatrixXd BD(m,maxvec);              // container for recursive product of B*D
  MatrixXd z_loc(maxvec,nsub);        //
  MatrixXd w_loc(m,nsub);             //

  // extra local containers
  ArrayXd  lamMax(nsub);   // container for computing max lambda
  ArrayXd  alphaMax(nsub); // container for computing max alpha
  ArrayXd  setMax(nsub);   // container for getting correct set
  ArrayXd  ymin(nsub);     // container for computing negative components of solution vector 

  Array<long int,Dynamic,1> jk(nsub); // container for getting correct index 

  // solution set containers
  std::vector<long int>               l(nsub);        // l[i] is the dimension of the subset of selected indices local to a subdomain.
  std::vector<std::vector<long int> > dualSet(nsub);  // local indices
  std::vector<std::vector<int> >      setKey(nsub);         // container for archiving which set an element belongs to (i.e. A or -A)
  std::list<std::pair<int,long_int> > gdualSet;       // global indices

  //initialization 
  rnorm = bnorm;

  for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);
 
  for(int i=0; i<nsub; ++i) {
    y[i] = VectorXd::Zero(maxlocvec[i]);
    g1[i].resize(A[i].cols());
    g2[i].resize(A[i].cols());
    g3[i].resize(A[i].cols());
    g4[i].resize(maxlocvec[i]);
    S[i].resize(A[i].cols());
    if(scaling) for(int j=0; j<A[i].cols(); ++j) S[i][j] = 1/A[i].col(j).norm();
    else S[i].setOnes();
  }
 
  r = b; 
  c.setZero();
  DtGDinv.setZero();
  alpha.setZero();
  lambda.setZero();

  for(int i=0; i<nsub; ++i) {
    B[i]  = MatrixXd::Zero(m,maxlocvec[i]);
    D[i]  = MatrixXd::Zero(maxlocvec[i],maxvec);
    GD[i] = MatrixXd::Zero(maxlocvec[i],maxvec);
  }
  
  // counters
  long int k      = 0; // cardinality of global solution set
  long int iter   = 0; // iteration counter
  long int downIt = 0; // downdate counter

/*
   char hostname[256];
   gethostname(hostname, sizeof(hostname));
   printf("PID %d on %s ready for attach\n", getpid(), hostname);
   fflush(stdout);
   sleep(25);
*/
 
  // ********************************* Main Loop ****************************************
  while (true) {

    if(myrank == 0 && verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "   "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k == maxvec || iter >= maxit) break;

    // compute step length along current face
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif 
    for(int i=0; i<nsub; ++i) {
      g1[i] = S[i].asDiagonal()*(A[i].transpose()*r);
      g2[i] = S[i].asDiagonal()*(A[i].transpose()*c);
      
      // make sure that element has not already been selected
       g3[i] = g1[i]; for(long int j=0; j<l[i]; ++j) g3[i][dualSet[i][j]] = 0.0;

      // compute maximum step length for this subdomain
      for(long int col = 0; col != A[i].cols(); col++){
        double num = g3[i][col];
        double den = g2[i][col];
        if(num >= 0.) {
          g3[i][col] = num/(1.0-den);
          g2[i][col] = std::numeric_limits<double>::min();
        } else {
          g3[i][col] = std::numeric_limits<double>::min();
          g2[i][col] = (-1.0)*num/(1.0+den);
        }
      }

      bool whichSet = 1;      
      long int position1 = 0;
      long int position2 = 0;

      double lam  = 0.0;
      double lam1 = g3[i].maxCoeff(&position1);
      double lam2 = g2[i].maxCoeff(&position2);

      if(A[i].cols() > 0){
        if (lam1 > lam2){
          lam = lam1;
          jk[i] = position1;
        } else {
          lam = lam2;
          jk[i] = position2;
          whichSet = -1;
        }
      } else {
        lam = 0.;
      }

      lamMax[i] = lam; 
      setMax[i] = whichSet;

    }

    int ik;  // subdomain which has max coeff
    int Set; // correct sign of selected set
    s.val   = lamMax.maxCoeff(&ik);
    s.rank  = myrank;
    p.index = jk[ik];
    p.sub   = ik; 
    Set     = setMax[ik];

#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
    MPI_Bcast(&Set, 1, MPI_INT, s.rank, mpicomm);
#endif

    lambda[k] = 1.0/s.val; 

    // add column to B and update GD, setkey and local dualSet
    if(s.rank == myrank) {
      B[ik].col(l[ik]) = S[ik][jk[ik]]*double(Set)*A[ik].col(jk[ik]);
      GD[ik].row(l[ik]).head(k) = B[ik].col(l[ik]).transpose()*BD.leftCols(k);
      dualSet[ik].push_back(jk[ik]);
      setKey[ik].push_back(Set);
      l[ik]++;
    }
    gdualSet.push_back(std::pair<int,long_int>(s.rank,p));

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      for(long int j=0; j<l[i]; ++j) g4[i][j] = g1[i][dualSet[i][j]];
      z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g4[i].head(l[i]);
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
      d.head(l[i]) = g4[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
      w_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
    }
    Block<MatrixXd,Dynamic,1,true> w = BD.col(k) = w_loc.rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, w.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
    DtGDinv[k] = 1/w.squaredNorm();
    alpha[k] = r.dot(w)*DtGDinv[k];
    r -= alpha[k]*w;
    c += lambda[k]*alpha[k]*w;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      y[i].head(l[i]) += alpha[k]*d.head(l[i]);
      GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*w;
    }
    k++;

    //****************************** Secondary Loop *******************************
    while (true) {
      iter++;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
      for(int i=0; i<nsub; ++i) {
        ymin[i] = (l[i] > 0) ? y[i].head(l[i]).minCoeff(&jk[i]) : std::numeric_limits<double>::max();
      }
      int ik;
      s.val  = ymin.minCoeff(&ik);
      s.rank = myrank;
#ifdef USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpicomm);
#endif

      if(s.val < 0) {
        downIt++;

        if(s.rank == myrank) {
          p.index = dualSet[ik][jk[ik]];
          p.sub = ik;
          std::vector<long int>::iterator pos1 = dualSet[ik].begin() + jk[ik];
          std::vector<int>::iterator      pos2 = setKey[ik].begin()  + jk[ik];
          dualSet[ik].erase(pos1);
          setKey[ik].erase(pos2);
          y[ik][l[ik]-1] = 0;
//          std::cout << "removing index " << p.index << " from subdomain " << ik << " on process with rank " << myrank << std::endl;
        }
 
#ifdef USE_MPI
        MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
        std::list<std::pair<int,long_int> >::iterator pos = std::find(gdualSet.begin(), gdualSet.end(), std::pair<int,long_int>(s.rank,p));
        std::list<std::pair<int,long_int> >::iterator fol = gdualSet.erase(pos);

        // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
        // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
        // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
        k = std::distance(gdualSet.begin(), fol);

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          y[i].head(l[i]).setZero();
          l[i] = 0; for(std::list<std::pair<int,long_int> >::iterator it = gdualSet.begin(); it != fol; ++it) { if(it->first == myrank && it->second.sub == i) l[i]++; }
          D[i].bottomRightCorner(maxlocvec[i]-l[i],maxvec-k).setZero(); // XXX this is a larger block than necessary
          y[i].head(l[i]) = D[i].topLeftCorner(l[i],k)*alpha.head(k);
        }
        r = b - BD.leftCols(k)*alpha.head(k); // XXX this could be parallelized, rowwise
        c = BD.leftCols(k)*lambda.head(k).asDiagonal()*alpha.head(k);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          g1[i].head(l[i]) = B[i].leftCols(l[i]).transpose()*r;
        }
 
        double num = 0.;
        double den = 0.;
        for(std::list<std::pair<int,long_int> >::iterator it = fol; it != gdualSet.end(); ++it) {
          if(it->first == myrank) {
            int i = it->second.sub;
            int slot;
            for(slot = 0; slot != dualSet[i].size(); slot++){
              if(it->second.index == dualSet[i][slot])
                break;
            }
            B[i].col(l[i]) = S[i][it->second.index]*double(setKey[i][slot])*A[i].col(it->second.index);
            g1[i][l[i]] = B[i].col(l[i]).dot(r);
            den = g1[i][l[i]];
            num = B[i].col(l[i]).dot(c);
            // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
            GD[i].row(l[i]).head(k) = B[i].col(l[i]).transpose()*BD.leftCols(k);
            l[i]++;
          }

#ifdef USE_MPI
          MPI_Bcast(&num, 1, MPI_DOUBLE, it->first, mpicomm);
          MPI_Bcast(&den, 1, MPI_DOUBLE, it->first, mpicomm);
#endif

          lambda[k] = (1.0-num)/den; 
//           std::cout << "num = " << num << " den = " << den << " rank = " << myrank << " lambda = " << lambda.head(k).transpose() << std::endl;

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g1[i].head(l[i]);
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
            d.head(l[i]) = g1[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
            w_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
          }
          Block<MatrixXd,Dynamic,1,true> w = BD.col(k) = w_loc.rowwise().sum();
#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, w.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

          DtGDinv[k] = 1/w.squaredNorm();
          alpha[k] = r.dot(w)*DtGDinv[k];
          r -= alpha[k]*w;
          c += lambda[k]*alpha[k]*w;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
            y[i].head(l[i]) += alpha[k]*d.head(l[i]);
            GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*w;
            g1[i].head(l[i]) -= alpha[k]*GD[i].col(k).head(l[i]);
          }
          k++;
        }
      } else { 
        break;
      } // end if block
    }

   rnorm = r.norm();

  }

  Array<VectorXd,Dynamic,1> x(nsub);
  for(int i=0; i<nsub; ++i) {
    x[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) x[i][dualSet[i][j]] += S[i][dualSet[i][j]]*double(setKey[i][j])*y[i][j];
  }
  return x;

}                   

#endif
