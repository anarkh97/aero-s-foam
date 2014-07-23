#include <cstdio>
#include <Math.d/Skyline.d/BlockSky.C>

template<>
void
GenBlockSky<complex<double> >::addImaginary(FullSquareMatrix &kel, int *dofs)
{
 // Construct stiffness matrix K (skyA)

 int i, j, ri, rj;

 int kndof = kel.dim();                 // Element stiffness dimension

 for( i = 0; i < kndof; ++i ) {          // Loop over rows.

   if( (ri = rowColNum[dofs[i]]) == -1 ) continue;// Skip constrained dofs

   for( j = 0; j < kndof; ++j ) {          // Loop over columns.

     if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.

     if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs

     skyA[dlp[rj] - rj + ri ] += complex<double>(0.0, kel[i][j]);
   }
 }

}

template<>
void
GenBlockSky<double>
   ::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenBlockSky<double> cannot addImaginary\n");
}

template<>
void
GenBlockSky<complex<double> >::add(FullSquareMatrixC &kel, int *dofs)
{
 // Construct stiffness matrix K (skyA)

 int i, j, ri, rj;

 int kndof = kel.dim();                 // Element stiffness dimension

 for( i = 0; i < kndof; ++i ) {          // Loop over rows.

   if( (ri = rowColNum[dofs[i]]) == -1 ) continue;// Skip constrained dofs

   for( j = 0; j < kndof; ++j ) {          // Loop over columns.

     if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.

     if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs

     skyA[dlp[rj] - rj + ri ] += kel[i][j];
   }
 }

}

template<>
void
GenBlockSky<double>
   ::add(FullSquareMatrixC &kel, int *dofs)
{
  fprintf(stderr, "GenBlockSky<double> cannot add FullSquareMatrixC\n");
}

template<>
void
GenBlockSky<double>::add(GenAssembledFullM<complex<double> > &kel, int *dofs)
{
 fprintf(stderr, "ERROR: add(GenAssembledFullM<complex<double> > &, int*) is being called"
                 "on an inappropriate matrix\n");
}

#define BLOCKSKY_INSTANTIATION_HELPER(Scalar) \
template \
GenBlockSky<Scalar>::~GenBlockSky(); \
template \
void \
GenBlockSky<Scalar>::initialize(); \
template \
GenBlockSky<Scalar>::GenBlockSky(Connectivity*, EqNumberer*, double); \
template \
GenBlockSky<Scalar>::GenBlockSky(Connectivity*, DofSetArray*, double); \
template \
GenBlockSky<Scalar>::GenBlockSky(Connectivity*, DofSetArray*, double, int*); \
template \
void \
GenBlockSky<Scalar>::factor(); \
template \
void \
GenBlockSky<Scalar>::parallelFactor(); \
template \
void \
GenBlockSky<Scalar>::add(FullSquareMatrix&, int*); \
template \
void \
GenBlockSky<Scalar>::add(int, int, Scalar); \
template \
void \
GenBlockSky<Scalar>::add(AssembledFullM&, int*); \
template \
void \
GenBlockSky<Scalar>::add(GenAssembledFullM<complex<double> >&, int*); \
template \
void \
GenBlockSky<Scalar>::add(FullM&, int, int); \
template \
void \
GenBlockSky<Scalar>::addBoeing(int, const int*, const int*, const double*, int*, Scalar); \
template \
void \
GenBlockSky<Scalar>::solve(GenVector<Scalar>&, GenVector<Scalar>&); \
template \
void \
GenBlockSky<Scalar>::solve(Scalar*, Scalar*); \
template \
void \
GenBlockSky<Scalar>::reSolve(Scalar*); \
template \
void \
GenBlockSky<Scalar>::reSolve(GenVector<Scalar>&); \
template \
void \
GenBlockSky<Scalar>::reSolve(int, Scalar**); \
template \
void \
GenBlockSky<Scalar>::print(FILE*); \
template \
void \
GenBlockSky<Scalar>::zeroAll(); \
template \
void \
GenBlockSky<Scalar>::clean_up(); \
template \
void \
GenBlockSky<Scalar>::unify(FSCommunicator*); \
template \
void \
GenBlockSky<Scalar>::addDiscreteMass(int, Scalar);

BLOCKSKY_INSTANTIATION_HELPER(double);
BLOCKSKY_INSTANTIATION_HELPER(complex<double>);

