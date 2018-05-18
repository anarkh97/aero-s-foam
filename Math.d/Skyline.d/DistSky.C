#include <Utils.d/dbg_alloca.h>
#include <cstdio>

#include <Comm.d/Communicator.h>
extern Communicator *structCom;

#include <Utils.d/linkfc.h>

#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, double *e, const int &f,
                  double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const complex<double> &d, complex<double> *e, const int &f,
                  complex<double> *g, const int &h, const complex<double> &i, complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

template<class Scalar>
GenDistSky<Scalar>::~GenDistSky()
{
  if(nlines) { delete [] nlines; nlines = 0; }
}

template<class Scalar>
GenDistSky<Scalar>::GenDistSky(Connectivity *cn, EqNumberer *dsa, 
                               double trbm, int fRow, int nRow) :
 GenSkyMatrix<Scalar>(cn, dsa, trbm)
{
 firstRow = fRow;
 numRows  = nRow;

 // allocate memory for the number of rows we wish to store
 nlines = new Scalar[numRows*this->numUncon];
}

template<class Scalar>
void
GenDistSky<Scalar>::parallelFactor()
{
 GenSkyMatrix<Scalar>::parallelFactor();

 // zero the n-lines
 int i;
 for(i=0; i<numRows*this->numUncon; ++i)
   nlines[i] = 0.0;   

 Scalar **rows = (Scalar **) dbg_alloca(sizeof(Scalar*)*numRows);

 for(i=0; i<numRows; ++i) {
   nlines[i*this->numUncon+firstRow+i] = 1.0;
   rows[i] = nlines+i*this->numUncon;
 }

 GenSkyMatrix<Scalar>::reSolve(numRows, rows);

 // delete entire Skyline matrix 
 delete [] this->skyA;
 this->skyA = 0;
}

template<class Scalar>
void
GenDistSky<Scalar>::reSolve(Scalar *rhs)
{
 Scalar *partialSum = (Scalar *)dbg_alloca(sizeof(Scalar)*numRows);

 Tgemv('T',this->numUncon,numRows,1.0,nlines,this->numUncon,
       rhs, 1, 0.0, partialSum, 1);

 // zero the rhs
 int i;
 for(i=0; i<this->numUncon; ++i)
   rhs[i] = 0.0;

 for(i=0; i<numRows; ++i)
   rhs[firstRow+i] = partialSum[i];

#ifdef DISTRIBUTED
  structCom->globalSum(this->numUncon, rhs);
#endif
}

template<class Scalar>
void
GenDistSky<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 if(this->neqs() > 0) reSolve(rhs.data());
}

template<class Scalar>
void
GenDistSky<Scalar>::reSolve(int nRHS, Scalar **RHS)
{
 if(this->neqs() > 0) GenSolver<Scalar>::reSolve(nRHS, RHS);
}

template<class Scalar>
void
GenDistSky<Scalar>::reSolve(int nRHS, GenVector<Scalar> *RHS)
{
 if(this->neqs() > 0) GenSolver<Scalar>::reSolve(nRHS, RHS);
}

template<class Scalar>
void
GenDistSky<Scalar>::reSolve(int nRHS, Scalar *RHS)
{
 if(this->neqs() > 0) GenSolver<Scalar>::reSolve(nRHS, RHS);
}

template<class Scalar>
void
GenDistSky<Scalar>::zeroAll()
{
 if(this->skyA == 0) this->skyA = new Scalar[GenSkyMatrix<Scalar>::size()];
 for(int i = 0; i < GenSkyMatrix<Scalar>::size(); ++i)
   this->skyA[i] = 0.0;
}

