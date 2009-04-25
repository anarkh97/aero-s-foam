#include <stdio.h>
#include <Math.d/Skyline.d/BlockSky.h>
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
