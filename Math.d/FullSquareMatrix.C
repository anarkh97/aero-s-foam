#include <iostream>
#include <stdio.h>
#include <cassert>
#include <algorithm>

#include <Math.d/FullSquareMatrix.h>

/*
extern "C"      {
   void _FORTRAN(dsyev)(const char &, const char &, const int &,
                        double *, const int &, double *, double *,
                        const int &, int &);
   //use zheev if complex matrix
   // for more functions see dsyevr and zheevr
}

inline void Tdsyev(const char &a, const char &b, const int &c,
                   double *d, const int &e, double *f, double *g,
                   const int &h, int &i)
{ _FORTRAN(dsyev)(a,b,c,d,e,f,g,h,i); }

*/
#include <Utils.d/MyComplex.h>

template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix()
{
 size = 0; value = 0; myval=1; length = 0;
}


template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix(int i, Scalar *l)
{
  size = i;
  length = size*size;
  if(l) {
    value=l;
    myval=0;
  }
  else {
    value=new Scalar[size*size];
    myval=1;
  }
}

template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix(GenFullSquareMatrix<Scalar> &m, Scalar s)
{
  size = m.dim();
  length = size*size;
  value= new Scalar[length];
  myval=1;
  int p = 0;
  for(int i=0; i<m.numRow(); ++i)
    for(int j=0; j<m.numCol(); ++j)
      value[p++] = m[i][j]*s;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::copy(const GenFullSquareMatrix<Scalar> &m)
{
  if(size != m.dim()) {
    size = m.dim();
    length = size*size;
    if(myval && value) delete [] value;
    value= new Scalar[length];
    myval=1;
  }
  int p = 0;
  for(int i=0; i<m.numRow(); ++i)
    for(int j=0; j<m.numCol(); ++j)
      value[p++] = m[i][j];
}

template<class Scalar>
GenFullSquareMatrix<Scalar>::~GenFullSquareMatrix()
{
  if(value && myval) delete [] value;
  value=0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::setSize(int i)
{
  if(value && myval) delete [] value;
  size = i;
  length = size*size;
  if(myval) value = new Scalar[length];
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::changeSize(int i, int numMax)
{
 if (i>numMax) fprintf(stderr," May be a Problem in GenFullSquareMatrix<Scalar>::changeSize ");
 else {
   size = i;
   length = size*size;
 }
}


template<class Scalar>
void
GenFullSquareMatrix<Scalar>::reSize(int newSize)
{
  if (newSize == size) {
    return;
  }
  
  // Create updated data
  int newLength = newSize * newSize;
  Scalar * newValue = new Scalar[newLength];
  int copySize = std::min(newSize, size);

  for (int i = 0; i < copySize; ++i) {
    const Scalar * rowBegin = value + i * size;
    Scalar * newRowBegin = newValue + i * newSize;
    std::copy(rowBegin, rowBegin + copySize, newRowBegin);
  } 

  // Commit changes
  size = newSize;
  length = newLength;
  if (myval) {
    delete[] value;
  }
  value = newValue;
  myval = 1;
}

template<class Scalar>
GenFullSquareMatrix<Scalar>
GenFullSquareMatrix<Scalar>::operator + (const GenFullSquareMatrix<Scalar> &M2)
{
  GenFullSquareMatrix<Scalar> res(size);
  for(int i = 0; i < length; ++i)
    res.setData(i, value[i] + M2.value[i]);
  return res;
}

template<class Scalar>
GenFullSquareMatrix<Scalar>
GenFullSquareMatrix<Scalar>::operator - (const GenFullSquareMatrix<Scalar> &M2)
{
  GenFullSquareMatrix<Scalar> res(size);
  for(int i = 0; i < length; ++i)
    res.setData(i, value[i] - M2.value[i]);
  return res;
}

template<class Scalar>
GenFullSquareMatrix<Scalar>
GenFullSquareMatrix<Scalar>::operator * (Scalar s)
{
  GenFullSquareMatrix<Scalar> res(size);
  for(int i = 0; i < length; ++i)
    res.setData(i, s*value[i]);
  return res;
}

template<class Scalar>
GenFullSquareMatrix<Scalar>
GenFullSquareMatrix<Scalar>::operator / (Scalar s)
{
  GenFullSquareMatrix<Scalar> res(size);
  for(int i = 0; i < length; ++i)
    res.setData(i, value[i]/s);
  return res;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator += (const GenFullSquareMatrix<Scalar> &M2)
{
  for(int i = 0; i < length; ++i)
    value[i] += M2.value[i];
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator -= (const GenFullSquareMatrix<Scalar> &M2)
{
  for(int i = 0; i < length; ++i)
    value[i] -= M2.value[i];
  return *this;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::multiply(GenFullSquareMatrix<Scalar> &M2, GenFullSquareMatrix<Scalar> &result)
{
  for(int i = 0; i < size; ++i) {
    for(int j = 0; j < size; ++j) {
      result[i][j] = 0.0;
      for(int k = 0; k < size; ++k) {
        result[i][j] += (*this)[i][k]*M2[k][j];
      }
    }
  }
}


template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator += (const Tensor_d2s0 &tens)
{
  for(int i = 0; i < length; ++i)
    value[i] += tens[i];
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator *= (Scalar s)
{
  for(int i = 0; i < length; ++i)
    value[i] *= s;
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator /= (Scalar s)
{
  for(int i = 0; i < length; ++i)
    value[i] /= s;
  return *this;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::zero()
{
  for(int i=0; i<length; ++i)
    value[i] = 0.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::unitary()
{
  for(int i=0; i<length; ++i)
    value[i] = 1.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::unitaryDiag()
{
  int i;
  for(i=0; i<length; ++i)
     value[i] = 0.0;
  for(i=0; i<size; ++i)
    (*this)[i][i] = 1.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::symmetrize()
{
 int i, j;
 for(i=0; i<size; ++i)
   for(j=0;  j<i; ++j)
     (*this)[j][i] = (*this)[i][j] = 0.5*((*this)[j][i] + (*this)[i][j]);
}

// ... THIS FUNCTION PRINTS A FULLSQUARE MATRIX
// ... IN MATHEMATICA MATRIX INPUT FORM
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::print(const char *msg, const char *msg2)
{
 int mathematica = 0;

 if(mathematica) {
   if(*msg) printf("%s\n",msg);
   int i,j;
   printf("{\n");
   for(i=0; i<size; ++i) {
     printf("{ ");
     for(j=0; j<size; ++j) {
       if(j==size-1)
         printf("%e ",(*this)[i][j]);
       else
         printf("%e, ",(*this)[i][j]);
     }
     if(i==size-1)
       printf("}");
     else
       printf("},");
     printf("\n");
   }
   printf("}\n");
 } else {
   if(*msg) fprintf(stderr,"%s\n",msg);
   int i,j;
   for(i=0; i<size; ++i) {
     for(j=0; j<size; ++j)
       //fprintf(stderr,"%s(%d,%d)=%3.2e, ",msg2,i+1,j+1,(*this)[i][j]);
       fprintf(stderr,"%3.2e,",(*this)[i][j]);
     fprintf(stderr,"\n");
   }
 }
}

template<>
inline
void
GenFullSquareMatrix<std::complex<double> >::print(const char*, const char*) {
  fprintf(stderr, "GenFullSquareMatrix<std::complex<double> >::print not implemented\n");
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::printDiagonals()
{
  int i;
  for(i=0; i<size; ++i)
    fprintf(stderr,"% e\n",(*this)[i][i]);
}


template<class Scalar>
void
GenFullSquareMatrix<Scalar>::copy(Scalar *d)
{
  // WARNING : No check is performed on the size of d
  for(int i=0; i<length; ++i)
    value[i] = d[i];
}
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::add(GenFullSquareMatrix<Scalar> &m, int *rc)
{
  for(int i=0; i<m.dim(); ++i)
    for(int j=0; j<m.dim(); ++j)
      (*this)[rc[i]][rc[j]] += m[i][j];
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::multiply(GenFullSquareMatrix<Scalar> &res, double d)
{
  for(int i = 0; i < length; ++i)
    res.value[i] = value[i] * d;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::eigenVals(Scalar* eigenV)
{
  //Scalar work[6*size];
  Scalar *copyValue = new Scalar[size*size];
  for (int i = 0 ; i < size*size ; i++) copyValue[i] = value[i];
  //use the lapack routine dsyev
  char a = 'N';//'N' for eigenVlaues alone and 'V' for eigenvalues + eigenvectors
  char b = 'U';//store upper triangle
  int c = size;//order of the matrix
  double *d = copyValue;
  int e = size;
  double *f = eigenV;
  double *g = new Scalar[6*size];
  int h = 6*size;
  int i = 0;
//  _FORTRAN(dsyev)(a, b, c, d, e, f, g, h, i);
  //other function can be used: zheev for complex or dsyevr & zheevr for faster methods
  Tdsyev(a, b, c, d, e, f, g, h, i);
  delete[] g;
  delete[] copyValue;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::eigenV(Scalar* eigenV)
{
  // value will change to store the eigenvectors
  //use the lapack routine dsyev
  char a = 'V';//'N' for eigenVlaues alone and 'V' for eigenvalues + eigenvectors
  char b = 'U';//store upper triangle
  int c = size;//order of the matrix
  double *d = value;
  int e = size;
  double *f = eigenV;
  double *g = new Scalar[6*size];
  int h = 6*size;
  int i = 0;
//  _FORTRAN(dsyev)(a, b, c, d, e, f, g, h, i);
  //other function can be used: zheev for complex or dsyevr & zheevr for faster methods
  Tdsyev(a, b, c, d, e, f, g, h, i);
  delete[] g;
}

/*
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::invert(GenFullSquareMatrix<Scalar> GFSM)
{
//see sommerElement.C
}
*/

#include <Utils.d/linkfc.h>
// BLAS prototypes
extern "C" {
  // BLAS level three real Matrix-VectorProduct
  extern void _FORTRAN(dgemv)(const char &, const int &,const int &,
			      const double &, double *, const int &,
			      double *, const int &, const double &, double *, const int &);

  extern void _FORTRAN(zgemv)(const char &, const int &,const int &,
			      const DComplex &, DComplex *, const int &,
			      DComplex *, const int &, const DComplex &, DComplex *, const int &);
}

template<> template<>
inline void GenFullSquareMatrix<DComplex>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, double c, GenFullSquareMatrix<DComplex>::TransposeFlag transposed)
{
  DComplex alpha = c;
  DComplex beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(zgemv)(trans, size, size, alpha, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<DComplex>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, DComplex c, GenFullSquareMatrix<DComplex>::TransposeFlag transposed)
{
  DComplex beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(zgemv)(trans, size, size, c, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<double>& a, GenVector<double>& b, double c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  double beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(dgemv)(trans, size, size, c, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, double c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  double beta  = 1.0;
  int two = 2;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(dgemv)(trans, size, size, c, value, size,
		  reinterpret_cast<double*>(a.data()), two, beta,
		  reinterpret_cast<double*>(b.data()), two);
  _FORTRAN(dgemv)(trans, size, size, c, value, size,
		  reinterpret_cast<double*>(a.data())+1, two, beta,
		  reinterpret_cast<double*>(b.data())+1, two);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<double>& a, GenVector<double>& b, DComplex c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  assert(0);
}

template<class Scalar> template<class Scalar1, class Scalar2, class Scalar3>
void GenFullSquareMatrix<Scalar>::multiply(GenVector<Scalar1>& a, GenVector<Scalar2>& b, Scalar3 c, typename GenFullSquareMatrix<Scalar>::TransposeFlag transposed)
{
  // multiply  b += c * matrix * a
  // Warning: no check for correct dimensions
  Scalar1* ap = a.data();
  Scalar2* bp = b.data();
  if (transposed == NORMAL) {
    for(int i=0;i<size;i++)
    {
      for(int j=0;j<size;j++)
      {
        bp[i] += c*(value+size*i)[j]*ap[j];
      }
    }
  } else {
    for(int i=0;i<size;i++)
    {
      for(int j=0;j<size;j++)
      {
        bp[j] += c*(value+size*i)[j]*ap[i];
      }
    }
  }
  return;
}
