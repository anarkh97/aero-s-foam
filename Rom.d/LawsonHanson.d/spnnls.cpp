#include <cmath>
#include <iostream>
#include <iomanip>
#ifdef USE_EIGEN3
#include <Eigen/Core>

//     SUBROUTINE SPNNLS  (A,MDA,M,N,B,X,RELTOL,RNORM,W,ZZ,ZZ2,INDEX,MODE)

//  Algorithm SPNNLS: SPARSE NONNEGATIVE LEAST SQUARES

//  The original version of this code was developed by
//  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
//  1973 JUN 15, and published in the book
//  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
//  Revised FEB 1995 to accompany reprinting of the book by SIAM.

//  Adapted by Julien Cortial at Stanford University

//     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE A
//     SPARSE N-VECTOR, X, THAT VERIFIES

//         ||A * X - B||_2 <= RELTOL * ||B||_2  SUBJECT TO X .GE. 0
//     ------------------------------------------------------------------
//                     Subroutine Arguments

//     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE
//                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N
//                     MATRIX, A.           ON EXIT A() CONTAINS
//                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
//                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY
//                     THIS SUBROUTINE.
//     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON-
//             TAINS Q*B.
//     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL
//             CONTAIN THE SOLUTION VECTOR.
//     RELTOL  RELATIVE TOLERANCE
//             (STOPPING CRITERION: ||B - A*X||_2 < RELTOL * ||B||_2).
//     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE
//             RESIDUAL VECTOR.
//     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN
//             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.
//             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z
//     ZZ()     AN M-ARRAY OF WORKING SPACE.
//     ZZ2()    AN N-ARRAY OF WORKING SPACE.
//     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
//                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
//                 P AND Z AS FOLLOWS..

//                 INDEX(1)   THRU INDEX(NSETP) = SET P.
//                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.
//                 IZ1 = NSETP + 1 = NPP1
//                 IZ2 = N
//     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING
//             MEANINGS.
//             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
//             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.
//                   EITHER M .LE. 0 OR N .LE. 0.
//             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS.

//     ------------------------------------------------------------------

template<typename AnyMatrix>
int spnnls(AnyMatrix &A, int mda, int m, int 
           n, double *_b, double *_x, double reltol, double 
           &rnorm, double *_w, double *_zz, double *_zz2, int *_index,
           int &mode)
{
    // Builtin functions
    using std::sqrt;
    using std::pow;
    using std::abs;
    using std::cout;
    using std::endl;
    using std::setw;
    using std::setprecision;

    // Local variables
    static int i__, j, l;
    static double t;
    extern int g1(double, double, double &, double &, double &);
    static double cc;
    extern int h12(int, int, int, int, double *, int, double &, double *, int, int, int);
    static int ii, jj, ip;
    static double sm;
    static int iz;
    static double up, ss;
    static int iz1, iz2, npp1;
    extern double diff(double, double);
    static int iter;
    static double temp, wmax, alpha, asave;
    static int itmax, izmax, nsetp;
    static double dummy, unorm, ztest, abstol;

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > b(_b, m);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > x(_x, n);

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > w(_w, n);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > zz(_zz, m);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > zz2(_zz2, n);
    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 1> > index(_index, n);

    Eigen::Matrix<double, Eigen::Dynamic, 1> a_j(m), a_jj(m);

    // Function Body
    mode = 1;
    if (m <= 0 || n <= 0) {
      mode = 2;
      return 0;
    }
    iter = 0;
    itmax = n * 3;

//  INITIALIZE THE ARRAYS INDEX AND X.
    x.setZero();
    index.setLinSpaced(n,0,n); // index[i] = i

    iz2 = n-1;
    iz1 = 0;
    nsetp = 0;
    npp1 = 1;

//  INIT ZZ2 = tr(A^T * A)^{-1}
    for (int i = 0; i < n; ++i) {
      zz2[i] = 0.;
      for (int l = 0; l < m; ++l) {
        zz2[i] += pow(A.coeff(l, i), 2);
      }
      zz2[i] = (double)1. / sqrt(zz2[i]);
    }

//  INIT ABSTOL = RELTOL^2 * ||B||^2
    sm = b.squaredNorm();
    abstol = reltol * sqrt(sm);

//  ******  MAIN LOOP BEGINS HERE  ******
    do {

//      QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
//      OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.
        if (iz1 > iz2 || nsetp >= m) {
            break;
        }
    
//      COMPUTE THE NORM^2 OF THE RESIDUAL VECTOR.
        sm = b.segment(npp1-1,m-npp1-1).squaredNorm();
        rnorm = sqrt(sm);
    
        cout << " Iteration= " << setw(11) << iter; 
        cout << " Active set size = " << setw(11) << nsetp;
        cout << " Residual norm = " << setw(20) << setprecision(16) << rnorm << "     ";
        cout << " Target = " << setw(20) << setprecision(16) << abstol << "     " << endl;
    
//      STOPPING CRITERION
        if (rnorm < abstol) {
          break;
        }
    
//      COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
        for (iz = iz1; iz < iz2+1; ++iz) {
          j = index[iz];
          sm = 0.;
          for (int l = npp1-1; l < m; ++l) {
            sm += A.coeff(l, j) * b[l];
          }
          w[j] = sm * zz2[j];
        }

        bool stop = false;
        do {
//        FIND LARGEST POSITIVE W(J).
          wmax = 0.;
          for (iz = iz1; iz < iz2+1; ++iz) {
            j = index[iz];
            if (w[j] > wmax) {
              wmax = w[j];
              izmax = iz;
            }
          }
    
//        IF WMAX .LE. 0. GO TO TERMINATION.
//        THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
          if (wmax <= 0.) {
            stop = true;
            break;
          }
          iz = izmax;
          j = index[iz];
    
//        THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
//        BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
//        NEAR LINEAR DEPENDENCE.
          asave = A.coeff(npp1-1, j);

          for(int row=0; row<m; ++row) a_j[row] = A.coeff(row,j);
          h12(1, npp1, npp1+1, m, a_j.data(), 1, up, &dummy, 1, 1, 0);
          for(int row=0; row<m; ++row) A.coeffRef(row,j) = a_j[row];

          unorm = 0.;
          if (nsetp != 0) {
            for (int l = 0; l < nsetp; ++l) {
              unorm += pow(A.coeff(l, j), 2);
            }
          }
          unorm = sqrt(unorm);
          if (diff(unorm + abs(A.coeff(npp1-1, j)) * .01, unorm) > 0.) {
    
//          COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
//          AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
            zz = b;

            for(int row=0; row<m; ++row) a_j[row] = A.coeff(row,j);
            h12(2, npp1, npp1+1, m, a_j.data(), 1, up, zz.data(), 1, 1, 1);

            ztest = zz[npp1-1] / A.coeff(npp1-1, j);
          } else  {
            ztest = 0;
          }
    
//        REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
//        RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
//        COEFFS AGAIN.
          if(ztest <= 0.) {
            A.coeffRef(npp1-1, j) = asave;
            w[j] = 0.;
          }
        } while(ztest <= 0.);
        if(stop) break;
    
//      THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
//      SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
//      TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
//      COL J,  SET W(J)=0.
        b = zz;
    
        index[iz] = index[iz1];
        index[iz1] = j;
        ++iz1;
        nsetp = npp1;
        ++npp1;
    
        if (iz1 <= iz2) {
          for (int jz = iz1; jz < iz2+1; ++jz) {
            jj = index[jz];

            for(int row=0; row<m; ++row) { a_j[row] = A.coeff(row,j); 
                                           a_jj[row] = A.coeff(row,jj); }
            h12(2, nsetp, npp1, m, a_j.data(), 1, up, a_jj.data(), 1, mda, 1);
            for(int row=0; row<m; ++row) A.coeffRef(row,jj) = a_jj[row];
          }
        }
    
        if (nsetp != m) {
          for (int l = npp1-1; l < m; ++l) {
            A.coeffRef(l, j) = 0.;
          }
        }
    
        w[j] = 0.;

//      SOLVE THE TRIANGULAR SYSTEM.
//      STORE THE SOLUTION TEMPORARILY IN ZZ().
        for (int l = 0; l < nsetp; ++l) {
          ip = nsetp - l;
          if (l != 0) {
            for (int ii = 0; ii < ip; ++ii) {
              zz[ii] -= A.coeff(ii, jj) * zz[ip];
            }
          }
          jj = index[ip-1];
          zz[ip-1] /= A.coeff(ip-1, jj);
        }
    
//      ******  SECONDARY LOOP BEGINS HERE ******
        ++iter;
        for( ; iter <= itmax; ++iter) {
        
//        SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
//        IF NOT COMPUTE ALPHA.
          alpha = 2.;
          for (ip = 0; ip < nsetp; ++ip) {
            int l = index[ip];
            if (zz[ip] <= 0.) {
              t = -x[l] / (zz[ip] - x[l]);
              if (alpha > t) {
                alpha = t;
                jj = ip;
              }
            }
          }
        
//        IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
//        STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.
          if (alpha == 2.) {
            break;
          }
        
//        OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
//        INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.
          for (ip = 0; ip < nsetp; ++ip) {
            int l = index[ip];
            x[l] += alpha * (zz[ip] - x[l]);
          }
        
//        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
//        FROM SET P TO SET Z.
          i__ = index[jj];
          do {
            x[i__] = 0.;
        
            if (jj != nsetp-1) {
              ++jj;
              for (int j = jj; j < nsetp; ++j) {
                ii = index[j];
                index[j - 1] = ii;
                g1(A.coeff(j - 1, ii), A.coeff(j, ii), cc, ss, A.coeffRef(j - 1, ii));
                A.coeffRef(j, ii) = 0.;
                for (int l = 0; l < n; ++l) {
                  if (l != ii) {

//                  Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
                    temp = A.coeff(j - 1, l);
               	    A.coeffRef(j - 1, l) = cc * temp + ss * A.coeff(j, l);
               	    A.coeffRef(j, l) = -ss * temp + cc * A.coeff(j, l);
               	  }
                }

//              Apply procedure G2 (CC,SS,B(J-1),B(J))
                temp = b[j - 1];
                b[j - 1] = cc * temp + ss * b[j];
                b[j] = -ss * temp + cc * b[j];
              }
            }
        
            npp1 = nsetp;
            --nsetp;
            --iz1;
            index[iz1] = i__;
 
//          SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
//          BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
//          IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
//          THAT ARE NONPOSITIVE WILL BE SET TO ZERO
//          AND MOVED FROM SET P TO SET Z.
            for (jj = 0; jj < nsetp; ++jj) {
              i__ = index[jj];
              if (x[i__] <= 0.) {
                break;
              }
            }
          } while (jj != nsetp);
        
//        COPY B INTO ZZ.  THEN SOLVE AGAIN AND LOOP BACK.
          zz = b;
          for (int l = 0; l < nsetp; ++l) {
            ip = nsetp - l;
            if (l != 0) {
              for (int ii = 0; ii < ip; ++ii) {
                zz[ii] -= A.coeff(ii, jj) * zz[ip];
              }
            }
            jj = index[ip-1];
            zz[ip-1] /= A.coeff(ip-1, jj);
          }
        }
//      ******  END OF SECONDARY LOOP  ******
    
        if (iter >= itmax) {
          mode = 3;
          break;
        }
    
        for (ip = 0; ip < nsetp; ++ip) {
          i__ = index[ip];
          x[i__] = zz[ip];
        }

//      ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.

    } while(true);
//  ******  END OF MAIN LOOP  ******

//  COME TO HERE FOR TERMINATION.
//  COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.
    sm = 0.;
    if (npp1 <= m) {
      for (int i = npp1-1; i < m; ++i) {
        sm += pow(b[i], 2);
      }
    } else {
      w.setZero();
    }
    rnorm = sqrt(sm);
    return 0;
}
#endif
