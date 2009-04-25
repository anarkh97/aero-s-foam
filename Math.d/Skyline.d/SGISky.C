#include <stdio.h>
//#include <fstream.h>
#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>

#include <Utils.d/linkfc.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Math.d/Skyline.d/utility.h>
#include <Corotational.d/GeomState.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/GetTime.h>

#ifdef sgi

extern "C" {
void _FORTRAN(dskyj)(double *values, int* pointers, double *dinv,
                     const int& ncpu, int& neq, const int& lead);
void _FORTRAN(fbsj)(double *values, double *dinv, int& n,
                    int* pointers, double *b, int *scratch);
}

#endif

SGISky::~SGISky()
{
  if(skyA) { delete skyA; skyA = 0; }
  if(dinv) { delete dinv; dinv = 0; }
}

void
SGISky::zeroAll()
{
  int i;
  for(i = 0; i < dlp[numUncon-1]; ++i)
    skyA[i] = 0.0;
}

void
SGISky::clean_up()
{
 if(skyA) {
   delete [] skyA;
   skyA = 0;
 }
}

void
SGISky::allMult(double x)
{
 int i;
 for(i = 0; i < dlp[numUncon-1]; ++i)
   skyA[i] *= x;
}

SGISky::SGISky(Connectivity *cn, EqNumberer *, ConstrainedDSA *c_dsa, 
                     double trbm, Rbm *rigid) :
SkyData(cn,c_dsa,trbm,rigid)
{
  constructTime = -getTime();

  // ALLOCATE MEMORY FOR SKYA 
  skyA = new double[ dlp[numUncon-1] ];

  // initialize to zero
  zeroAll();

  dinv = new double[dlp[numUncon-1]];

  solveTime = 0.0;

  constructTime += getTime();

}

SGISky::SGISky(Connectivity *cn, EqNumberer *_dsa, double trbm, int *bc) :
SkyData(cn,_dsa,trbm,bc)
{
  constructTime = -getTime();

  // ALLOCATE MEMORY FOR SKYA 
  skyA = new double[ dlp[numUncon-1] ];

  // INITIALIZE SKYA TO ZERO
  zeroAll();

  dinv = new double[dlp[numUncon-1]];

  solveTime = 0.0;

  constructTime += getTime();

}

SGISky::SGISky(Connectivity *cn, EqNumberer *_dsa, double trbm, int *rCN, int) :
SkyData(_dsa,cn,trbm,rCN)
{
  constructTime = -getTime();

  // ALLOCATE MEMORY FOR SKYA   
  skyA = new double[ dlp[numUncon-1] ];

  // INITIALIZE SKYA TO ZERO
  zeroAll();

  dinv = new double[dlp[numUncon-1]];

  solveTime = 0.0;

  constructTime += getTime();

}


SGISky::SGISky(FullM *mat, double tolerance) :
SkyData(mat->numRow(), tolerance)
{

  constructTime = -getTime();

  // WARNING: DEFAULT value of tolerance is 1.0E-4

  if (mat->numRow() != mat->numCol()) {
    fprintf(stderr, " ERROR: SGISky::SGISky(FullM *mat) constructor\n");
    skyA = 0;
  }
  else {
    skyA = new double[ dlp[numUncon-1] ];

    int offset = 0;
    int i,j;
    for(i = 0; i < numUncon; ++i)
      for(j = 0; j <= i; ++j) {
        skyA[offset++] = (*mat)[i][j];
      }
  }

  dinv = new double[dlp[numUncon-1]];

  constructTime += getTime();
  
}

void
SGISky::printMemory()
{
   fprintf(stderr,"Memory necessary for Skyline array (Mb): %10.3f\n",
          8*dlp[numUncon-1]/(1024.0*1024.0));
}

void
SGISky::printConstructTime()
{
  fprintf(stderr,"SGISky Construction Time = %14.5f\n",constructTime/1000.0);
}

void 
SGISky::mult(const Vector &, Vector &)
{
  fprintf(stderr,"This shouldn't be called--SGISky::mult\n");
}

void
SGISky::factor()
{
 #ifdef sgi
  _FORTRAN(dskyj)(skyA, dlp, dinv, 6, numUncon, 1);
 #else
   fprintf(stderr,"SGI SKYLINE ONLY AVAILABLE ON SGI\n");
 #endif
}

void
SGISky::solve(Vector &rhs, Vector &solution)
{
   solution = rhs;

   reSolve(solution.data());
}

void
SGISky::solve(double *rhs, double *solution)
{
   int i;
   for(i=0; i<dim(); ++i)
     solution[i] = rhs[i];

   reSolve(solution);
}


void
SGISky::reSolve( Vector &rhs )
{
  reSolve(rhs.data()); // ML Need to uniformize resolve
}

void
SGISky::reSolve( double *rhs )
{

 solveTime -= getTime();


 #ifdef sgi
   int * scr = (int *) dbg_alloca(sizeof(int)*numUncon);
   _FORTRAN(fbsj)(skyA, dinv, numUncon, dlp, rhs, scr);
 #else
   fprintf(stderr,"SGI SKYLINE ONLY AVAILABLE ON SGI\n");
 #endif

 solveTime += getTime();

}

void
SGISky::reSolve(int nRHS, double **rhs)
{
 int i;
 for(i=0; i<nRHS; ++i)
   reSolve(rhs[i]);
}

void
SGISky::reSolve(int nRHS, Vector *rhs)
{
 int i;
 for(i=0; i<nRHS; ++i)
   reSolve(rhs[i].data());
}

void
SGISky::getNullSpace(double *rbm)
{
}

Vector *
SGISky::getNullSpace()
{
 return 0;
}

void
SGISky::getRBMs(double *rigidBodyModes)
{
 rbm->getRBMs(rigidBodyModes);
}

void
SGISky::getRBMs(Vector *rigidBodyModes)
{
 rbm->getRBMs(rigidBodyModes);
}

void
SGISky::getRBMs(VectorSet& rigidBodyModes)
{
 rbm->getRBMs(rigidBodyModes);
}

double
SGISky::diag(int dof) const
{
  if(skyA[ dlp[dof] - 1]  == 0 )
    return (1.0);
  else
    return skyA[ dlp[dof] - 1];
}

double &
SGISky::diag(int dof)
{
  return skyA[ dlp[dof] - 1];
}

void
SGISky::add(FullSquareMatrix &kel, int *dofs)
{
// Construct stiffness matrix K (skyA)
 int i, j, k;
 int kndof = kel.dim();                	// Element stiffness dimension
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( rowColNum[dofs[i]] == -1 ) continue;    // Skip constrained dofs 
    for( j = 0; j < kndof; ++j ) {       	// Loop over columns.
       if( dofs[i] > dofs[j] ) continue;        // Work with upper symmetric half.
       if( rowColNum[dofs[j]] == -1 ) continue; // Skip constrained dofs
       k = dlp[rowColNum[dofs[j]]] - (rowColNum[dofs[j]] - rowColNum[dofs[i]]) - 1;
       skyA[k] += kel[i][j];
    }
 }
}

void
SGISky::add(FullM &knd, int fRow, int fCol)
{
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

  int iCol, iRow;
  for(iCol = 0; iCol < ncol; ++iCol) {
    int fk = dlp[fCol+iCol] + fRow - fCol - iCol - 1;
    int rowStop = (nrow < fCol+iCol-fRow+1) ? nrow : fCol+iCol-fRow+1;
    for(iRow = 0; iRow < rowStop; ++iRow)
      skyA[fk+iRow] += knd[iRow][iCol];
  }
}

/*
void 
SGISky::addDiscreteMass(int cdof, double dmass)
{
  skyA[dlp[cdof]-1] += dmass;
}
*/

void
SGISky::addDiscreteMass(int dof, double dmass)
{
  if(dof < 0) return;
  int cdof = rowColNum[dof];
  if(cdof < 0) return;
  skyA[dlp[cdof]-1] += dmass;
}


