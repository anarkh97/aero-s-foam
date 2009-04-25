#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Utils.d/Connectivity.h>
#include <Utils.d/linkfc.h>
#include <Math.d/Vector.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/SGISparseMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/GetTime.h>

#ifdef sgi
extern "C" {
 void _FORTRAN(psldlt_preprocess)(int &, int &, int *, int *, int &, double *);
 void _FORTRAN(psldlt_factor)(int &, int &, int *, int *, double *);
 void _FORTRAN(psldlt_solve)(int &, double *, double *);

 void _FORTRAN(zpsldlt_preprocess)(int &, int &, int *, int *, int &, DComplex *);
 void _FORTRAN(zpsldlt_factor)(int &, int &, int *, int *, DComplex *);
 void _FORTRAN(zpsldlt_solve)(int &, DComplex *, DComplex *);
}

inline void Tpsldlt_preprocess(int &a, int &b, int *c, int *d, int &e, double *f)
 { _FORTRAN(psldlt_preprocess)(a,b,c,d,e,f); }
inline void Tpsldlt_preprocess(int &a, int &b, int *c, int *d, int &e, DComplex *f)
 { _FORTRAN(zpsldlt_preprocess)(a,b,c,d,e,f); }

inline void Tpsldlt_factor(int &a, int &b, int *c, int *d, double *e)
 { _FORTRAN(psldlt_factor)(a,b,c,d,e); }
inline void Tpsldlt_factor(int &a, int &b, int *c, int *d, DComplex *e)
 { _FORTRAN(zpsldlt_factor)(a,b,c,d,e); }

inline void Tpsldlt_solve(int &a, double *b, double *c)
 { _FORTRAN(psldlt_solve)(a,b,c); }
inline void Tpsldlt_solve(int &a, DComplex *b, DComplex *c)
 { _FORTRAN(zpsldlt_solve)(a,b,c); }
#endif

template<class Scalar>
GenSGISparseMatrix<Scalar>::GenSGISparseMatrix(Connectivity *cn, DofSetArray *_dsa, 
                                               ConstrainedDSA *c_dsa, Rbm *_rbm, int _token) :
  SparseData(_dsa,c_dsa,cn)
{
  token = _token;

  rbm = _rbm;
  if(rbm) 
    numrbm = rbm->numRBM();
  else
    numrbm = 0;

  // Allocate memory for SGI sparse matrix
  unonz = new Scalar[xunonz[numUncon]];

  zeroAll();

#ifdef sgi
  // Time symbolic factorization done in psldlt_preprocess
  double constructTime = - getTime();

//  fprintf(stderr,"SGI Sparse Matrix Memory: \n");
//  fprintf(stderr,"(12 bytes X ) %25lld %14.3f Mb\n",
//          size(),12.0*size()/(1024.0*1024.0));
//  fprintf(stderr,"( 4 bytes X ) %25d %14.3f Mb\n",
//          numUncon,4.0*numUncon/(1024.0*1024.0));

  Tpsldlt_preprocess(token, numUncon, xunonz, rowu, nonz, &ops);

//  fprintf(stderr,"( 8 bytes X ) %25d %14.3f Mb\n",
//                  nonz,8.0*nonz/(1024.0*1024.0));

  constructTime += getTime();
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

  newrhs = 0;
  solveTime = 0.0;
}

template<class Scalar>
GenSGISparseMatrix<Scalar>::~GenSGISparseMatrix()
{
 if(unonz) { delete [] unonz; unonz = 0; }
 if(newrhs) { delete [] newrhs; newrhs = 0; }
}

template<class Scalar>
double
GenSGISparseMatrix<Scalar>::getMemoryUsed()
{
 double  essential_memory_size =
 //(xunonz[numUncon] * 12.0 +  nonz * 8.0 +  numUncon * 4.0) / (1024.0*1024.0);
 ((long int)xunonz[numUncon] * 12.0 + (long int)numUncon*4.0) / (1024.0*1024.0);

 return essential_memory_size;
}

template<class Scalar>
long
GenSGISparseMatrix<Scalar>::size()
{
 return (long) xunonz[numUncon];
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::factor()
{
#ifdef sgi
  Tpsldlt_factor(token, numUncon, xunonz, rowu, unonz);
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::solve(Scalar *rhs, Scalar *solution)
{
 solveTime -= getTime();

#ifdef sgi
 Tpsldlt_solve(token, solution, rhs);
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

 solveTime += getTime();
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution)
{
 solve(rhs.data(), solution.data());
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 solveTime -= getTime();

#ifdef sgi
 Tpsldlt_solve(token, rhs.data(), rhs.data());
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

 solveTime += getTime();
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::reSolve(Scalar *rhs)
{
 solveTime -= getTime();

#ifdef sgi
 Tpsldlt_solve(token, rhs, rhs);
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

 solveTime += getTime();
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::reSolve(int numRHS, Scalar **RHS)
{
 solveTime -= getTime();

#ifdef sgi
 int i;
 for(i=0; i<numRHS; ++i)
   Tpsldlt_solve(token, RHS[i], RHS[i]);
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

 solveTime += getTime();
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::reSolve(int numRHS, GenVector<Scalar> *RHS)
{
 solveTime -= getTime();

#ifdef sgi
 int i;
 for(i=0; i<numRHS; ++i)
   Tpsldlt_solve(token, RHS[i].data(), RHS[i].data());
#else
  fprintf(stderr, "SGI SPARSE SOLVER IS SUPPORTED ONLY ON SGI MACHINES\n");
  exit(-1);
#endif

 solveTime += getTime();
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::zeroAll()
{
  int i;
  for(i=0; i < xunonz[numUncon]; ++i)
    unonz[i] = 0.0;
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::clean_up()
{
 if(unonz) {
   delete [] unonz;
   unonz = 0;
 }
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::add(FullSquareMatrix &kel, int *dofs)
{
 int i, j, m, mstart, mstop;

 int kndof = kel.dim();                          // Dimension of element stiffness.
 for(i = 0; i < kndof; ++i ) {                   // Loop over rows.
    if(unconstrNum[dofs[i]] == -1 ) continue;    // Skip constrained dofs
    for(j = 0; j < kndof; ++j ) {                // Loop over columns.
       if(dofs[i] > dofs[j] ) continue;          // Work with upper symmetric half.
       if(unconstrNum[dofs[j]] == -1 ) continue; // Skip constrained dofs
       mstart = xunonz[unconstrNum[dofs[j]]];
       mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m = mstart; m < mstop; ++m) {
          if(rowu[m-1] == (unconstrNum[dofs[i]] + 1) ) {
            unonz[m-1] += kel[i][j];
            break;
          }
       }
       // if(m == mstop) fprintf(stderr," *** ERROR: SGISparseMatrix::add \n" );
    }
 }
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::print(const char *message)
{
 if(message) fprintf(stderr,"%s\n",message);
 int i;
 for(i=0; i < xunonz[numUncon]; ++i)
    fprintf(stderr,"%e\n",unonz[i]);
}

template<class Scalar>
Scalar
GenSGISparseMatrix<Scalar>::diag(int dof) const
{
  int mstart = xunonz[dof]-1;
  int mstop  = xunonz[dof+1]-1;

  int m;
  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
      if(unonz[m] == 0.0)
        return (1.0);
      else
        return unonz[m];
    }
  }
  throw "GenSGISparseMatrix<Scalar>::diag - 1 - this should never be reached";
}

template<class Scalar>
Scalar &
GenSGISparseMatrix<Scalar>::diag(int dof)
{
  int mstart = xunonz[dof]-1;
  int mstop  = xunonz[dof+1]-1;

  int m;
  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
        return unonz[m];
    }
  }
  throw "GenSGISparseMatrix<Scalar>::diag - 2 - this should never be reached";
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
 if(dof < 0) return;
 int cdof = unconstrNum[dof];
 if(cdof < 0) return;
 int diagLocator = xunonz[dof+1]-2; // This should be the diagonal   
 unonz[diagLocator] += dmass;
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::getRBMs(double *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::getRBMs(Vector *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSGISparseMatrix<Scalar>::getRBMs(VectorSet& rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

