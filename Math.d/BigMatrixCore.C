#include <Utils.d/linkfc.h>
#include <Math.d/BigMatrix.h>

// BLAS level three real Matrix Product

extern "C" {

  void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                       const int &, const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemm)(const char &, const char &, const int &,const int &,
                       const int &, const complex<double> &, complex<double> *, const int &,
                       complex<double> *, const int &, const complex<double> &, complex<double> *,                       const int &);

  void _FORTRAN(dgemv)(const char &, const int &,const int &,
                       const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemv)(const char &, const int &,const int &,
                       const complex<double> &, complex<double> *, const int &,
                       complex<double> *, const int &, const complex<double> &, complex<double> *, const int &);
  void _FORTRAN(dsisl)(double *, const int &, const int &, int *, double *);
  void _FORTRAN(dsifa)(double *, const int &, const int &, int *, int &);
}


template<>
void
GenBigMatrix<double>::lineUpdate(int iThread, int numThreads, int iBlock)
{
 int rem = iBlock%numThreads;
 int myBlock = iBlock-rem+iThread;
 if(myBlock <= iBlock) myBlock += numThreads;
 int i;

 for(;myBlock < numBlocks; myBlock += numThreads) {
   // symmetric copy
   symCopy(block(myBlock,iBlock), block(iBlock,myBlock),
        blockSize(myBlock),blockSize(iBlock));
   // multiply Rji^T by Rii^-1
   for(i = 0; i < blockSize(myBlock); ++i)
   _FORTRAN(dsisl)(block(iBlock,iBlock), size, blockSize(iBlock),
        kpvt[iBlock],block(myBlock,iBlock)+i*size);
 }
}


template<>
void
GenBigMatrix<DComplex>::lineUpdate(int iThread, int numThreads, int iBlock)
{
 fprintf(stderr,"GenBigMatrix<DComplex>::lineUpdate(int iThread, int numThreads, int iBlock) not implemented\n");
}

template<>
void
GenBigMatrix<double>::factorDiagonal(int iThread, int numThreads, int iBlock)
{
 int info;
 if(iThread == 0)
  _FORTRAN(dsifa)(block(iBlock,iBlock),size,blockSize(iBlock), kpvt[iBlock],info);
}

template<>
void
GenBigMatrix<DComplex>::factorDiagonal(int iThread, int numThreads, int iBlock)
{
 fprintf(stderr, "GenBigMatrix<DComplex>::factorDiagonal(int iThread, int numThreads, int iBlock) not implemented\n");
}

template<>
void
GenBigMatrix<double>::diagSolve(int iBlock, int iThread, int numThreads, double *rhs)
{
 if(iThread == 0)
 _FORTRAN(dsisl)(block(iBlock,iBlock), size, blockSize(iBlock),
        kpvt[iBlock],rhs+blockIndex[iBlock]);
}

template<>
void
GenBigMatrix<DComplex>::diagSolve(int iBlock, int iThread, int numThreads, DComplex *rhs)
{
 fprintf(stderr, "GenBigMatrix<DComplex>::diagSolve(int iBlock, int iThread, int numThreads, DComplex *rhs) not implemented \n");
}

