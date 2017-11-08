#include <Math.d/BLKSparseMatrix.C>

template<> 
void
GenBLKSparseMatrix<complex<double> >
   ::addImaginary(FullSquareMatrix &kel, int *dofs)
{
 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();

 for(i = 0; i < kndof; ++i ) {
   if((rowi = unconstrNum[dofs[i]]) == -1) continue;
   p1 = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper = invsuper[p1] - 1;
   fstcol = xsuper[csuper] - 1;
   lxbeg = xlindx[csuper]-1;
   lxend = xlindx[csuper+1]-1;
   for(j = 0; j < kndof; ++j ) {
     if((colj = unconstrNum[dofs[j]]) == -1) continue;
     int irow   = invp[colj] - 1;
     if(irow >= fstcol) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[position - 2 - offset] += complex<double>(0.0, kel[i][j]);
           break;
         }
       }
     }
   }
 }
}

template<> 
void
GenBLKSparseMatrix<complex<double> >
   ::add(FullSquareMatrixC &kel, int *dofs)
{

 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();

 for( i = 0; i < kndof; ++i ) {
   if( (rowi = unconstrNum[dofs[i]]) == -1 ) continue;
     p1     = invp[rowi] - 1;
     position = xlnz[p1+1];
     csuper = invsuper[p1] - 1;
     fstcol = xsuper[csuper] - 1;
     lxbeg  = xlindx[csuper]-1;
     lxend  = xlindx[csuper+1]-1;
     for( j = 0; j < kndof; ++j ) {
       if( (colj = unconstrNum[dofs[j]]) == -1 ) continue;
       int irow   = invp[colj] - 1;
       if ( irow >= fstcol ) {
         offset = lxend - lxbeg;
         for(k=lxbeg; k<lxend; ++k) {
           offset -= 1;
           if(lindx[k]-1 == irow) {
             lnz[ position - 2 - offset] += kel[i][j];
             break;
           }
         }
       }

     }
 }

}

template<> 
void
GenBLKSparseMatrix<double>
   ::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenBLKSparseMatrix<double> cannot addImaginary\n");
}

template<> 
void
GenBLKSparseMatrix<double>
   ::add(FullSquareMatrixC &kel, int *dofs)
{
  fprintf(stderr, "GenBLKSparseMatrix<double> cannot add FullSquareMatrixC\n");
}

#define BLKSPARSEMATRIX_INSTANTIATION_HELPER(Scalar) \
template \
GenBLKSparseMatrix<Scalar>::~GenBLKSparseMatrix(); \
template \
void \
GenBLKSparseMatrix<Scalar>::init(); \
template \
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(Connectivity*, DofSetArray*, DofSetArray*, double, SolverCntl&, Rbm*); \
template \
void \
GenBLKSparseMatrix<Scalar>::factor(); \
template \
void \
GenBLKSparseMatrix<Scalar>::computeRBMs(); \
template \
void \
GenBLKSparseMatrix<Scalar>::getNullSpace(Scalar*); \
template \
double \
GenBLKSparseMatrix<Scalar>::getMemoryUsed(); \
template \
void \
GenBLKSparseMatrix<Scalar>::solve(Scalar*, Scalar*); \
template \
void \
GenBLKSparseMatrix<Scalar>::solve(GenVector<Scalar>&, GenVector<Scalar>&); \
template \
void  \
GenBLKSparseMatrix<Scalar>::reSolve(Scalar*); \
template \
void \
GenBLKSparseMatrix<Scalar>::reSolve(GenVector<Scalar>&); \
template \
void \
GenBLKSparseMatrix<Scalar>::reSolve(int, Scalar**); \
template \
void \
GenBLKSparseMatrix<Scalar>::reSolve(int, GenVector<Scalar>*); \
template \
void \
GenBLKSparseMatrix<Scalar>::zeroAll(); \
template \
void \
GenBLKSparseMatrix<Scalar>::clean_up(); \
template \
void \
GenBLKSparseMatrix<Scalar>::unify(FSCommunicator*); \
template \
void \
GenBLKSparseMatrix<Scalar>::add(FullSquareMatrix&, int*); \
template \
void \
GenBLKSparseMatrix<Scalar>::add(GenAssembledFullM<Scalar>&, int*); \
template \
void \
GenBLKSparseMatrix<Scalar>::add(FullM&, int, int); \
template \
void \
GenBLKSparseMatrix<Scalar>::add(Scalar*); \
template \
void \
GenBLKSparseMatrix<Scalar>::addBoeing(int, const int*, const int*, const double*, int*, Scalar); \
template \
void \
GenBLKSparseMatrix<Scalar>::addone(Scalar, int, int); \
template \
void \
GenBLKSparseMatrix<Scalar>::mult(const Scalar*, Scalar*); \
template \
Scalar \
GenBLKSparseMatrix<Scalar>::getone(int, int); \
template \
void \
GenBLKSparseMatrix<Scalar>::print(); \
template \
void \
GenBLKSparseMatrix<Scalar>::printAll(); \
template \
Scalar \
GenBLKSparseMatrix<Scalar>::diag(int) const; \
template \
Scalar & \
GenBLKSparseMatrix<Scalar>::diag(int); \
template \
void \
GenBLKSparseMatrix<Scalar>::addDiscreteMass(int, Scalar); \
template \
void \
GenBLKSparseMatrix<Scalar>::getRBMs(double*); \
template \
void \
GenBLKSparseMatrix<Scalar>::getRBMs(Vector*); \
template \
void \
GenBLKSparseMatrix<Scalar>::getRBMs(VectorSet&); \
template \
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(Connectivity*, DofSetArray*, int*, double, SolverCntl&); \
template \
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(Connectivity*, EqNumberer*, double, SolverCntl&, int); \
template \
void \
GenBLKSparseMatrix<Scalar>::allocateMemory(); \
template \
int \
GenBLKSparseMatrix<Scalar>::numRBM(); \
template \
void \
GenBLKSparseMatrix<Scalar>::reSolve(GenFullM<Scalar>*);

BLKSPARSEMATRIX_INSTANTIATION_HELPER(double);
BLKSPARSEMATRIX_INSTANTIATION_HELPER(complex<double>);

