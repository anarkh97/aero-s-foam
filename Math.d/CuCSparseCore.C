#include <Math.d/CuCSparse.h>
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

