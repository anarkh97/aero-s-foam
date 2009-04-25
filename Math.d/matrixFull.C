#include <Math.d/matrix.h>
template<>
void
GenAssembledFullM<complex<double> >::addImaginary(FullSquareMatrix &mat, int *dofs)
{
 int kndof = mat.dim();
 int i,j,ri,cj;
 for(i=0; i<kndof; ++i) {
   if(dofs[i] == -1) continue;
   if((ri = rowMap[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(dofs[j] == -1) continue;
     if( (cj = colMap[dofs[j]]) == -1) continue;
     (*this)[ri][cj] += complex<double>(0.0, mat[i][j]);
   }
 }
}

template<>
void
GenAssembledFullM<double>::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenAssembledFullM<double> cannot addImaginary\n");
}

template<>
void
GenAssembledFullM<complex<double> >::add(FullSquareMatrixC &mat, int *dofs)
{
 int kndof = mat.dim();
 int i,j,ri,cj;
 for(i=0; i<kndof; ++i) {
   if(dofs[i] == -1) continue;
   if((ri = rowMap[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(dofs[j] == -1) continue;
     if( (cj = colMap[dofs[j]]) == -1) continue;
     (*this)[ri][cj] += mat[i][j];
   }
 }
}

template<>
void
GenAssembledFullM<double>::add(FullSquareMatrixC &kel, int *dofs)
{
  fprintf(stderr, "GenAssembledFullM<double> cannot add FullSquareMatrixC\n");
}

