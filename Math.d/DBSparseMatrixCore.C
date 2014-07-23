#include <cstdio>
#include <Math.d/DBSparseMatrix.C>

template<>
void
GenDBSparseMatrix<complex<double> >::add(FullSquareMatrixC &kel, int *dofs)
{
 int i, j, m;
 int kndof = kel.dim();                 // Dimension of element stiffness.
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( unconstrNum[dofs[i]] == -1 ) continue;      // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {              // Loop over columns.
       if( dofs[i] > dofs[j] ) continue; // Work with upper symmetric half.
       if( unconstrNum[dofs[j]] == -1 ) continue;   // Skip constrained dofs
       int mstart = xunonz[unconstrNum[dofs[j]]];
       int mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m=mstart; m<mstop; ++m)
       {
          if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) )
          {
            unonz[m-1] += kel[i][j];
            break;
          }
       }
       // if(m == mstop) fprintf(stderr," *** ERROR: GenDBSparseMatrix<Scalar>::add\n");
    }
 }
}

template<>
void
GenDBSparseMatrix<double>::add(FullSquareMatrixC &kel, int *dofs)
{
 fprintf(stderr, "Error: adding complex to real matrix\n");
}

template<>
void
GenDBSparseMatrix<complex<double> >::addImaginary(FullSquareMatrix &kel, int *dofs)
{
 int i, j, m;
 int kndof = kel.dim();                 // Dimension of element stiffness.
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( unconstrNum[dofs[i]] == -1 ) continue;      // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {              // Loop over columns.
       if( dofs[i] > dofs[j] ) continue; // Work with upper symmetric half.
       if( unconstrNum[dofs[j]] == -1 ) continue;   // Skip constrained dofs
       int mstart = xunonz[unconstrNum[dofs[j]]];
       int mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m=mstart; m<mstop; ++m)
       {
          if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) )
          {
            unonz[m-1] += complex<double>(0.0, kel[i][j]);
            break;
          }
       }
       // if(m == mstop) fprintf(stderr," *** ERROR: GenDBSparseMatrix<Scalar>::add\n");
    }
 }
}

template<>
void
GenDBSparseMatrix<double>
   ::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenDBSparseMatrix<double> cannot addImaginary\n");
}

/*
template<>
void
GenDBSparseMatrix<double>::transposeMult(const Vector &rhs, Vector &result)
{
  _FORTRAN(sptmv)(unonz, xunonz, rowu, numUncon, rhs.data(), result.data() );
}

template<>
void
GenDBSparseMatrix<DComplex>::transposeMult(const ComplexVector &rhs, ComplexVector &result)
{
  fprintf(stderr, "GenDBSparseMatrix<DComplex>::transposeMult(const ComplexVector &rhs, ComplexVector &result) not implemented \n");
}
*/

template<>
void
GenDBSparseMatrix<double>::multcomplex(const DComplex *rhs, DComplex *result)
{
 int nn = numUncon;
  _FORTRAN(cdspsmvp)(nn, unonz, xunonz, rowu, rhs, result );
}

template<>
void
GenDBSparseMatrix<DComplex>::multcomplex(const DComplex *rhs, DComplex *result)
{
  mult(rhs, result);
}

#define DBSPARSE_INSTANTIATION_HELPER(Scalar) \
template \
GenDBSparseMatrix<Scalar>::~GenDBSparseMatrix(); \
template  \
void \
GenDBSparseMatrix<Scalar>::clean_up(); \
template  \
void \
GenDBSparseMatrix<Scalar>::zeroAll(); \
template  \
void \
GenDBSparseMatrix<Scalar>::print(); \
template \
void \
GenDBSparseMatrix<Scalar>::print(char*); \
template  \
void \
GenDBSparseMatrix<Scalar>::print1(int, FILE*); \
template \
GenFullM<Scalar> * \
GenDBSparseMatrix<Scalar>::getFullMatrix(); \
template \
void \
GenDBSparseMatrix<Scalar>::add(int, int, Scalar); \
template  \
void \
GenDBSparseMatrix<Scalar>::add(FullSquareMatrix&, int*); \
template  \
void \
GenDBSparseMatrix<Scalar>::add(GenFullM<Scalar>&, int, int); \
template  \
void \
GenDBSparseMatrix<Scalar>::addBoeing(int, const int*, const int*, const double*, int*, Scalar); \
template  \
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(Connectivity*, DofSetArray*, int*); \
template  \
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(Connectivity*, DofSetArray*, ConstrainedDSA*); \
template  \
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(Connectivity*, EqNumberer*); \
template  \
double \
GenDBSparseMatrix<Scalar>::getMemoryUsed(); \
template  \
void \
GenDBSparseMatrix<Scalar>::mult(const Scalar*, Scalar*); \
template \
void \
GenDBSparseMatrix<Scalar>::mult(const GenVector<Scalar>&, Scalar*);  \
template \
void \
GenDBSparseMatrix<Scalar>::transposeMult(const GenVector<Scalar>&, GenVector<Scalar>&); \
template \
void \
GenDBSparseMatrix<Scalar>::transposeMult(const Scalar*, Scalar*); \
template \
void \
GenDBSparseMatrix<Scalar>::multAdd(const Scalar*, Scalar*); \
template  \
void \
GenDBSparseMatrix<Scalar>::mult(const GenVector<Scalar>&, GenVector<Scalar>&); \
template  \
Scalar \
GenDBSparseMatrix<Scalar>::diag(int) const; \
template  \
Scalar & \
GenDBSparseMatrix<Scalar>::diag(int); \
template  \
void \
GenDBSparseMatrix<Scalar>::makeIdentity(); \
template  \
void \
GenDBSparseMatrix<Scalar>::invertDiag(); \
template \
void \
GenDBSparseMatrix<Scalar>::addDiscreteMass(int, Scalar); \
template  \
long \
GenDBSparseMatrix<Scalar>::size(); \
template  \
void \
GenDBSparseMatrix<Scalar>::multDiag(const Scalar*, Scalar*); \
template  \
void \
GenDBSparseMatrix<Scalar>::multDiag(int, const Scalar**, Scalar**); \
template  \
void \
GenDBSparseMatrix<Scalar>::unify(FSCommunicator*); \
template \
void \
GenDBSparseMatrix<Scalar>::symmetricScaling(); \
template \
void \
GenDBSparseMatrix<Scalar>::applyScaling(Scalar*);

DBSPARSE_INSTANTIATION_HELPER(double);
DBSPARSE_INSTANTIATION_HELPER(complex<double>);
