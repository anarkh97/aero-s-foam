#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <algorithm>
#include <stdio.h>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>
#include <set>

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
// A       = where A*x= b
// b       = target vector 
// rnorm   = residual norm 
// maxsze  = maximum allowable sparsity
// reltol  = stopping criteria
// scaling = flag to turn on/off unit normalization of columns of A
//
Eigen::VectorXd
gpfp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      double maxsze, double reltol, bool verbose, bool scaling)
{
  using namespace Eigen;
  
  const long int columns = A.cols();				     // number of atoms in dictionary
  const long int m       = b.rows();                                 // number of measurements 
  const long int maxvec  = std::min(m, (long int)(maxsze*A.cols())); // maximum problem size
  const long int maxit   = 3*A.cols();                               // max allowable iterations
  double bnorm           = b.norm();        			     // norm of target
  double abstol          = reltol*bnorm;			     // stopping criteria

  VectorXd direction(maxvec);  // working array for  A
  VectorXd gradient(maxvec);   // working array for -A
  VectorXd y(maxvec);          // working array
  VectorXd r(m);               // residual
  VectorXd c(m);               // vertex estimate
  VectorXd update(m);          // working array
  VectorXd g1(columns);        // gradient
  VectorXd g2(columns);        // gradient working array 
  VectorXd DtGDinv(maxvec);    // inverse of D^T*A^T*A*D (see ref. 2)
  VectorXd lambda(maxvec);     // step lengths for vertex estimate updates
  VectorXd alpha(maxvec);      // step lengths for solution and residual updates 
  VectorXd S(columns);         // scaling to normalize columns of A
 
  MatrixXd B(m,maxvec);        // storage for selected columns of [A,-A]
  MatrixXd D(maxvec,maxvec);   // storage for conjugate directions
  MatrixXd GD(maxvec,maxvec);  // storage for D^T*G (see blumensath)

  std::vector<long int> dualSet;    // container for indices of dual set of [A,-A]
  std::vector<int>      setKey;     // container for archiving which set an element belongs to

  //initialization
  direction.setZero();  
  gradient.setZero();
  y.setZero();

  c.setZero();
  r     = b;
  rnorm = bnorm;
  update.setZero();

  B.setZero(); 
  D.setZero(); 
  GD.setZero();

  long int k      = 0; // number of selected elements
  long int iter   = 0;        // number of iterations
  long int downIt = 0;        // number of downdates

  // vestigial scaling 
  if(scaling) for(int i=0; i<A.cols(); ++i) S[i] = 1/A.col(i).norm();
  else S.setOnes();

  // *********************** Main Loop ********************************
  while (true) {
    // Output Iteration info
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
    bool set1 = true;
    long int position  = 0;
    long int position1 = 0;   
    long int position2 = 0;

    // compute step length along current face
    g1 = S.asDiagonal()*(A.transpose()*r); 
    g2 = S.asDiagonal()*(A.transpose()*c);

    //zero out elements that are already selected
    for(long int j=0; j<k; ++j) g1[dualSet[j]] = 0.0;

    for(long int col = 0; col != columns; col++){
      double num = g1(col); 
      double den = g2(col);
      if(num > 0) {
        g1(col) = num/(1-den);
        g2(col) = 0.;
      } else {
        g1(col) = 0.;
        g2(col) = (-1)*num/(1+den);
      }
    }

    double lam  = 0.0;
    double lam1 = g1.maxCoeff(&position1); 
    double lam2 = g2.maxCoeff(&position2);

    if (lam1 > lam2){
      lam = lam1;
      position = position1;
    } else {
      lam = lam2;
      position = position2;
      set1 = false;
    }

    long int curRowCol = k;
    long int numRowCol = k+1;

    // remove maximum element from the correct active set and place
    // in correct dual set, also store appropriate column of [A,-A]
    if(set1){ // element selected from A
      setKey.push_back(1);
      B.col(curRowCol) = S[position]*A.col(position);
    } else {  // element selected from -A
      setKey.push_back(-1);
      B.col(curRowCol) = S[position]*(-1.)*A.col(position);
    }
    dualSet.push_back(position); 

    // add step length for vertex estimate update to array
    lambda(k) = lam;

    // compute conjugate direction
    if (k == 0) {
      direction.head(numRowCol) = B.col(curRowCol).transpose()*r;
      GD.col(curRowCol).head(numRowCol) = direction.head(numRowCol).transpose()*B.col(curRowCol).transpose()*B.col(curRowCol);
    } else {
      // compute gradient vector
      gradient.head(numRowCol) = B.leftCols(numRowCol).transpose()*r;

      // update container for D^T*G, recall that it is upper triangular by construction
      GD.col(curRowCol).head(numRowCol-1) = D.topLeftCorner(curRowCol,curRowCol).triangularView<Upper>().transpose()*B.leftCols(numRowCol-1).transpose()*B.col(curRowCol);
 
      // now update the direction vector with one gigantic line....
      direction.head(numRowCol) = gradient.head(numRowCol) - 
          D.topLeftCorner(numRowCol,curRowCol).triangularView<Upper>()*(DtGDinv.head(curRowCol).asDiagonal()*(GD.topLeftCorner(curRowCol,numRowCol)*gradient.head(numRowCol)));
    }
 
    // update solution vector
    update  = B.leftCols(numRowCol)*direction.head(numRowCol);
    DtGDinv(curRowCol) = 1.0/(update.squaredNorm());

    alpha(curRowCol)   = r.dot(update)*DtGDinv(curRowCol);
    y.head(numRowCol) += alpha(curRowCol)*direction.head(numRowCol);

    // update residual
    r -= alpha(curRowCol)*update;
    
    // update vertex estimate
    c += alpha(curRowCol)*update/(lambda(curRowCol));

    D.col(curRowCol).head(numRowCol) = direction.head(numRowCol);

    rnorm = r.norm();
    
    // increment set size
    k++;

    // ****************************** Secondard Loop ************************************
    while (true) {
      iter++;
      // if any coefficients of y are negative, choose the minimum and remove it
      // all subsequent search direction must be re A-orthogonalized
      long int krestart = 0;
      double minVal = y.head(numRowCol).minCoeff(&krestart);
      if (minVal < 0) {
      downIt++;
  
/*    FOR DEBUGGING
      char hostname[256];
      gethostname(hostname, sizeof(hostname));
      printf("PID %d on %s ready for attach\n", getpid(), hostname);
      fflush(stdout);
      sleep(30);
*/

      // update element set
      dualSet.erase(dualSet.begin()+krestart);
      setKey.erase(setKey.begin()+krestart);

      // decrement set size
      k--;

      // rebuild residual and vertex estimate
      r = b;
      c.setZero();
      for(int i = 0; i < krestart; i++){
        g1 = B.leftCols(i+1)*D.col(i).head(i+1);
        c += alpha(i)*g1/(lambda(i));       
        r -= alpha(i)*g1;
      }

      //re A-orthogonalize search directions starting from krestart
      //all previous direction are already A-orthogonal
      for(int j = krestart; j < k+1; j++){

        int curRowCol = j;
        int numRowCol = j+1;
 
        // update container for selected columns of A
        B.col(curRowCol) = S[dualSet[j]]*setKey[j]*A.col(dualSet[j]);

        direction.setZero();
        if (curRowCol == 0) {
          direction.head(numRowCol) = B.col(curRowCol).transpose()*r;
          GD.col(curRowCol).head(numRowCol) = direction.head(numRowCol).transpose()*B.col(curRowCol).transpose()*B.col(curRowCol);
        } else {
          // compute gradient vector
          gradient.head(numRowCol) = B.leftCols(numRowCol).transpose()*r;
      
          // update container for D^T*G, recall that it is upper triangular by construction
          GD.col(curRowCol).head(numRowCol-1) = D.topLeftCorner(curRowCol,curRowCol).triangularView<Upper>().transpose()*B.leftCols(numRowCol-1).transpose()*B.col(curRowCol);
      
          // now update the direction vector with one gigantic line....
          direction.head(numRowCol) = gradient.head(numRowCol) -
                      D.topLeftCorner(numRowCol,curRowCol).triangularView<Upper>()*(DtGDinv.head(curRowCol).asDiagonal()*(GD.topLeftCorner(curRowCol,numRowCol)*gradient.head(numRowCol)));
        } 
       
        // update solution vector
        update  = B.leftCols(numRowCol)*direction.head(numRowCol);
        DtGDinv(curRowCol) = 1.0/(update.squaredNorm());
 
        // update step directions       
        alpha(curRowCol)   = r.dot(update)*DtGDinv(curRowCol);
        lambda(curRowCol)  = (B.col(curRowCol).dot(r))/(1-B.col(curRowCol).dot(c));

 
        // update residual
        r -= alpha(curRowCol)*update;
        
        // update vertex estimate
        c += alpha(curRowCol)*update/(lambda(curRowCol));
       
        // update direction container 
        D.col(curRowCol).head(numRowCol) = direction.head(numRowCol);
  
       }//end re A-orthogonlization
    
      y.setZero(); 
      y.head(k+1) +=  D.topLeftCorner(k+1,k+1)*alpha.head(k+1);
      for(int q = 0; q != m; q++) if(abs(y(q)) < 1e-6) y(q) = 0.0;
      rnorm = r.norm();     

      } else {
        break;
      }

    } // end secondary loop

  } // end main loop

  if(verbose) std::cout.flush();

  // reconstruct the solution
  VectorXd x = VectorXd::Zero(columns);
  int element  = 0;
  for(std::vector<int>::iterator it = setKey.begin(); it != setKey.end(); it++) {
    x[dualSet[element]] +=  S[dualSet[element]]*(*it)*y[element];
    element++;
  }  
  return x;
}

#endif
