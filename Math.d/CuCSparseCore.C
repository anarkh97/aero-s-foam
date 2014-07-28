#include <Math.d/CuCSparse.C>

template<>
void
GenCuCSparse<complex<double> >::add(FullSquareMatrixC &kel, int *dofs)
{
 int i, j, m, ri;
 int kndof = kel.dim();
 for(i=0; i<kndof; ++i) {
   if((ri = unconstrNum[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(constrndNum[dofs[j]] == -1) continue;
     for(m=xunonz[constrndNum[dofs[j]]]; m<xunonz[constrndNum[dofs[j]]+1]; ++m) {
       if(rowu[m] == ri) {
         Kuc[m] += kel[i][j];
         break;
       }
     }
   }
 }
}

template<>
void
GenCuCSparse<double>::add(FullSquareMatrixC &kel, int *dofs)
{
  fprintf(stderr, "GenCuCSparse<double> cannot add FullSquareMatrixC \n");
}

template<>
void
GenCuCSparse<complex<double> >::addImaginary(FullSquareMatrix &kel, int *dofs)
{
 int i, j, m, ri;
 int kndof = kel.dim();

 for(i=0; i<kndof; ++i) {
   if((ri = unconstrNum[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(constrndNum[dofs[j]] == -1) continue;
     for(m=xunonz[constrndNum[dofs[j]]]; m<xunonz[constrndNum[dofs[j]]+1]; ++m) {
       if(rowu[m] == ri) {
         Kuc[m] += complex<double>(0.0, kel[i][j]);
         break;
       }
     }
   }
 }
}

template<>
void
GenCuCSparse<double>::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenCuCSparse<double> cannot addImaginary\n");
}

#define CUCSPARSE_INSTANTIATION_HELPER(Scalar) \
template \
GenCuCSparse<Scalar>::~GenCuCSparse();  \
template  \
void \
GenCuCSparse<Scalar>::clean_up(); \
template  \
void \
GenCuCSparse<Scalar>::zeroAll(); \
template  \
void \
GenCuCSparse<Scalar>::negate(); \
template  \
void \
GenCuCSparse<Scalar>::print(); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(Connectivity*, DofSetArray*, int*); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(Connectivity*, DofSetArray*, DofSetArray*); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(Connectivity*, DofSetArray*, int*, int*); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(LMPCons**, int, DofSetArray*); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(int, int*, int, Scalar*, int); \
template  \
GenCuCSparse<Scalar>::GenCuCSparse(int, int, int*, int*, Scalar*); \
template  \
void  \
GenCuCSparse<Scalar>::doWeighting(Scalar*);  \
template  \
void  \
GenCuCSparse<Scalar>::doWeighting(int*); \
template  \
double \
GenCuCSparse<Scalar>::getMemoryUsed(); \
template  \
void \
GenCuCSparse<Scalar>::multSubtract(const GenVector<Scalar>&, GenVector<Scalar>&); \
template  \
void \
GenCuCSparse<Scalar>::multSubtract(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::mult(const GenVector<Scalar>&, GenVector<Scalar>&); \
template  \
void \
GenCuCSparse<Scalar>::mult(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::transposeMult(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::add(FullSquareMatrix&, int*); \
template  \
void \
GenCuCSparse<Scalar>::multIdentity(Scalar**); \
template  \
void \
GenCuCSparse<Scalar>::multIdentity(Scalar*); \
template \
void \
GenCuCSparse<Scalar>::multIdentity(Scalar**, int, int);  \
template  \
void \
GenCuCSparse<Scalar>::multSub(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::transposeMultSubtract(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::transposeMultSubtractClaw(const Scalar*, Scalar*, int, int*); \
template  \
void \
GenCuCSparse<Scalar>::multSub(int, Scalar**, Scalar**); \
template  \
void \
GenCuCSparse<Scalar>::multAdd(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::transposeMultAdd(const Scalar*, Scalar*); \
template  \
void \
GenCuCSparse<Scalar>::addBoeing(int, const int*, const int*, const double*, int*, Scalar); \
template  \
long \
GenCuCSparse<Scalar>::size(); \
template \
void \
GenCuCSparse<Scalar>::addDiscreteMass(int, Scalar); \
template \
void \
GenCuCSparse<Scalar>::add(int, int, Scalar); \
template \
void \
GenCuCSparse<Scalar>::multSubAdd(Scalar*, Scalar*); \
template \
void \
GenCuCSparse<Scalar>::multAddNew(const Scalar*, Scalar*); \
template \
void \
GenCuCSparse<Scalar>::transposeMultNew(const Scalar*, Scalar*); \
template \
void \
GenCuCSparse<Scalar>::transposeMultAddNew(const Scalar*, Scalar*); \
template \
void \
GenCuCSparse<Scalar>::transposeMultSubNew(const Scalar*, Scalar*); \
template \
void \
GenCuCSparse<Scalar>::mult(const Scalar*, Scalar*, Scalar, Scalar); \
template \
void \
GenCuCSparse<Scalar>::trMult(const Scalar*, Scalar*, Scalar, Scalar);

CUCSPARSE_INSTANTIATION_HELPER(double);
CUCSPARSE_INSTANTIATION_HELPER(complex<double>);

