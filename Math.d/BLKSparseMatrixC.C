#include <stdio.h>
#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>

#include <Utils.d/Memory.h>
#include <Utils.d/linkfc.h>
#include <Element.d/Element.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/matrixC.h>
#include <Math.d/BLKSparseMatrixC.h>
#include <Timers.d/GetTime.h>

extern "C" {

void _FORTRAN(ordmmd2)(int& n, int* xadj, int* adj, int* invp, int* perm,
                      int& iwsize, int* iwork, int& nofsub, int& iflag);

void _FORTRAN(sfinit)(int& n,   int& nnza,  int* xadj,  int* adj,   int* perm,
                      int* invp,int& maxsup,int& defblk,int* colcnt,int& nnzl,
                      int& nsub,int& nsuper,int* xsuper,int* snode, int& iwsize,
                      int* iwork, int& iflag);

void _FORTRAN(symfct)(int& n, int& nnza, int* xadj, int* adj, int* perm,
                      int* invp, int* colcnt, int& nsuper, int* xsuper, 
                      int* snode, int& nofsub, int* xlindx, int* lindx, 
                      int* xlnz, int& iwsize, int* iwork, int& iflag);

void _FORTRAN(bfinit)(int& nsuper, int* xsuper, int* snode, int* xlindx, 
                      int* lindx, int& tmpsiz, int& rwsize);

void _FORTRAN(zblkldl)(int& nsuper, int* xsuper, int *pnode, int* xlindx,
                      int* lindx,  int* xlnz, DComplex *lnz, int& defblk,
                      int &asdef,  int& numZEM, int& lbdef, int *def,
                      double& tol,
                      int *iprow, int* ipcol, int& tmpsiz, DComplex *temp,
                      int& iwsize, int* iwork, int& rwsize, DComplex *rwork,
                      int& iflag);

void _FORTRAN(zblkslv)(int &nsuper, int* xsuper, int* xlindx, int *lindx,
                      int* xlnz,
                      DComplex *lnx, int& defblk, int& numZEM, int& lbdef,
                      int* def, int* iprow, int* ipcol, int* perm, int* invp,
                      DComplex *rhs, DComplex *sol, DComplex *temp);

void _FORTRAN(zblkslvp)(int &nsuper, int* xsuper, int* xlindx, int *lindx,
                      int* xlnz,
                      DComplex *lnx, int& defblk, int& numZEM, int& lbdef,
                      int* def, int* iprow, int* ipcol, int* perm, int* invp,
                      int &nrhs, DComplex *rhs, int &ldr, DComplex *sol, 
                      DComplex *temp);

void _FORTRAN(zblkslv2)(int &nsuper, int* xsuper, int* xlindx, int *lindx,
                      int* xlnz,
                      DComplex *lnx, int& defblk, int& numZEM, int& lbdef,
                      int* def, int* iprow, int* ipcol, int* perm, int* invp,
                      DComplex *r1, DComplex *r2,
                      DComplex *s1, DComplex *s2,
                      DComplex *t1, DComplex *t2);

void _FORTRAN(zblkslv3)(int &nsuper, int* xsuper, int* xlindx, int *lindx,
                      int* xlnz,
                      DComplex *lnx, int& defblk, int& numZEM, int& lbdef,
                      int* def, int* iprow, int* ipcol, int* perm, int* invp,
                      DComplex *r1, DComplex *r2, DComplex *r3,
                      DComplex *s1, DComplex *s2, DComplex *s3,
                      DComplex *t1, DComplex *t2, DComplex *t3);

void _FORTRAN(zblkslv4)(int &nsuper, int* xsuper, int* xlindx, int *lindx,
                      int* xlnz,
                      DComplex *lnx, int& defblk, int& numZEM, int& lbdef,
                      int* def, int* iprow, int* ipcol, int* perm, int* invp,
                      DComplex *r1, DComplex *r2, DComplex *r3, DComplex *r4,
                      DComplex *s1, DComplex *s2, DComplex *s3, DComplex *s4,
                      DComplex *t1, DComplex *t2, DComplex *t3, DComplex *t4);

void _FORTRAN(zblkns)(int &nsuper, int *xsuper, int *xlindx, int *lindx,
                     int *xlnz, DComplex *lnz, int &defblk, int &nrbm,
                     int &lbdef, int *def, int *ipcol, int *invp, DComplex *ns,
                     int &numUncon, DComplex *temp);

void _FORTRAN(zzeromat)(int &numUncon, int &nsuper, int *xsuper, 
                       int *xlnz, DComplex *lnz, int *invsuper);

}


BLKSparseMatrixC::BLKSparseMatrixC(Connectivity *cn, DofSetArray *_dsa, 
                                 ConstrainedDSA *c_dsa, double _tol) :
  SparseData(_dsa,c_dsa,cn,1) {

  tol   = _tol;
  numrbm = 0;

  allocateMemory();

}


void
BLKSparseMatrixC::factor() {

  if (iwork==0)
     iwork = new int[7*numUncon+3];

  int iflag; // Error flag 

//       ***************************************************
//       Numerical input into data structure for sparse LDL'
//       factorization.
//       ***************************************************
//
//       --------------------------------------------------------
//       INPNV ...   input numerical values into data structures.
//
//       Input:      N, COLPTR, ROWIDX, VALUES, PERM, INVP,
//                   NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
//       Output:     LNZ
//       Work:       IWORK(N)
//       --------------------------------------------------------

//       ************************
//       Numerical factorization.
//       ************************
//       ---------------------------------------------------
//       BFINIT ...  initialization for block factorization.
//
//       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX
//       Output:     TMPSIZ, RWSIZE
//       ---------------------------------------------------

  int rwsize;


  _FORTRAN(bfinit)(nsuper, xsuper, snode, xlindx,
                   lindx,    tmpsiz, rwsize );

  // TESTING:
  tmpsiz = tmpsiz * 2;
  //tmpsiz = tmpsiz * 4;

  DComplex *tmpvec = new DComplex[tmpsiz];
  //TESTING
  //rwsize*=2;
  DComplex *rwork  = new DComplex[rwsize];

  iwsiz  =  3 * numUncon + 2 * nsuper;

//       -------------------------------------------------------
//       ZBLKLDL ...  numerical factorization.
//
//       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX, XLNZ,
//                   LNZ, DEFBLK, TOL, TMPSIZ, IWSIZE, RWSIZE
//       Output:     LNZ, NDEF, LBDEF, DEF, IPROW, IPCOL, IFLAG
//       Work:       TMPVEC(TMPSIZ), IWORK(2*N+2*NSUPER),
//                   RWORK(RWSIZE)
//       -------------------------------------------------------


  // if numrbm > 0, lbdef = number of rbm in last block

  // perform full pivoting on last 10x10 block of sparse matrix
  // def = integer array, dimension numUncon, it identifies
  //       the columns of the matrix that are linearly dependent.

  int *deftemp = (int *) dbg_alloca(sizeof(int)*numUncon);

  iprow = 0;
  ipcol = 0;
  lbdef = 0;

  if (defblk > 0) {
    iprow = new int[defblk]; // contains row pivoting sequence for last block.
    ipcol = new int[defblk]; // contains col pivoting sequence for last block.
    lbdef = numrbm;
  }

  _FORTRAN(zblkldl)(nsuper, xsuper,  snode, xlindx,  lindx,
                    xlnz, lnz, defblk, lbdef, numrbm,  lbdef, 
                    deftemp,    tol,  iprow,  ipcol, tmpsiz, 
                    tmpvec, iwsiz,  iwork, rwsize,  rwork, iflag);

  if (iflag != 0)
    fprintf(stderr, "Error during sparse factor %d\n",iflag);

  if (numrbm != lbdef) {
    fprintf(stderr,"Num rbm = %d last block (size %d) %d tol %e\n",numrbm,
             defblk, lbdef,tol); 
    if (defblk>0) 
       numrbm = lbdef;
  }

//     IFLAG =  0: successful factorization.
//     IFLAG = 31: insufficient work space in tmpvec.
//     IFLAG = 32: insufficient work space in iwork.
//     IFLAG = 33: insufficient work space in rwork.

  def = 0;

  if(numrbm > 0) {
    def = new int[numrbm];
    int i;
    for(i=0; i<numrbm; ++i)
      def[i] = deftemp[i];
  }

  delete [] rwork;  rwork  = 0;
  delete [] tmpvec; tmpvec = 0;
  delete [] iwork;  iwork  = 0;

  dbg_alloca(0);

}


double
BLKSparseMatrixC::getMemoryUsed() {

 // Figure this out later!
 return 0.0;

}


void
BLKSparseMatrixC::solve(DComplex *rhs, DComplex *solution)
{

 solveTime -= getTime();

 DComplex *temp = (DComplex *) dbg_alloca(sizeof(DComplex)*numUncon);

//       ********************
//       Triangular solution.
//       ********************
//
//       --------------------------------------------------------------
//       ZBLKSLV ...  numerical forward/backward solution.
//
//       Input:      nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk,
//                   numrbm, lbdef, def, iprow, ipcol, perm, invp, solution
//       Output:     solution
//       Work:       temp[numUncon]
//       --------------------------------------------------------------

 _FORTRAN(zblkslv)(nsuper, xsuper, xlindx, lindx, xlnz, 
                     lnz, defblk, numrbm, lbdef,  def, 
                   iprow,  ipcol,   perm,  invp,  rhs, 
                solution,   temp);

 solveTime += getTime();

}


void
BLKSparseMatrixC::solve(ComplexVector &rhs, ComplexVector &solution )
{
 solve( rhs.data(), solution.data() );
}


void 
BLKSparseMatrixC::reSolve(DComplex *rhs)
{
 solveTime -= getTime();

 DComplex *temp     = (DComplex *) dbg_alloca(sizeof(DComplex)*numUncon);
 DComplex *solution = (DComplex *) dbg_alloca(sizeof(DComplex)*numUncon);

 _FORTRAN(zblkslv)(nsuper, xsuper, xlindx, lindx, xlnz, 
                     lnz, defblk, numrbm, lbdef,  def, 
                   iprow,  ipcol,   perm,  invp,  rhs, 
                solution,   temp);

 // copy solution back to the rhs vector
 int i;
 for (i=0; i < numUncon; i++)
   rhs[i] = solution[i];

 solveTime += getTime();
}


void
BLKSparseMatrixC::reSolve(ComplexVector &rhs) {

 reSolve(rhs.data());

}


void
BLKSparseMatrixC::reSolve(int nRHS, DComplex **RHS) {

 solveTime -= getTime();

 int i = 0;
 int multiple = 4;
 int j;

 if(nRHS <= 4)
   multiple = nRHS;

 DComplex *t1 = new DComplex[numUncon];
 DComplex *t2 = new DComplex[numUncon];
 DComplex *t3 = new DComplex[numUncon];
 DComplex *t4 = new DComplex[numUncon];

 DComplex *s1 = new DComplex[numUncon];
 DComplex *s2 = new DComplex[numUncon];
 DComplex *s3 = new DComplex[numUncon];
 DComplex *s4 = new DComplex[numUncon];

 switch (multiple) {
   default:
   case 4:
     {
       // zblkslv4 = forward/backward routine for 4 RHS vectors
       for( ; i < nRHS-3; i += 4) {
         _FORTRAN(zblkslv4)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   RHS[i+2], RHS[i+3], s1, s2, s3, s4, t1, t2, t3, t4);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
            RHS[i+2][j] = s3[j];
            RHS[i+3][j] = s4[j];
         }
       }
     }
   case 3:
     {
       // zblkslv3 = forward/backward routine for 3 RHS vectors
       for( ; i < nRHS-2; i += 3) {
         _FORTRAN(zblkslv3)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   RHS[i+2], s1, s2, s3, t1, t2, t3);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
            RHS[i+2][j] = s3[j];
         }
       }
     }
   case 2:
     {
       // zblkslv2 = forward/backward routine for 2 RHS vectors
       for( ; i < nRHS-1; i += 2) {
         _FORTRAN(zblkslv2)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   s1, s2, t1, t2);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
         }
       }
     }
   case 1:
     {
       // zblkslv = forward/backward routine for 1 RHS vectors
       for ( ; i < nRHS; ++i) {
         _FORTRAN(zblkslv)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i],
                   s1,   t1); 
         for (j=0; j<numUncon; j++)
            RHS[i][j] = s1[j];
       }
     }
     break;
 }

 delete [] t1; t1 = 0;
 delete [] t2; t2 = 0;
 delete [] t3; t3 = 0;
 delete [] t4; t4 = 0;

 delete [] s1; s1 = 0;
 delete [] s2; s2 = 0;
 delete [] s3; s3 = 0;
 delete [] s4; s4 = 0;

 solveTime += getTime();

}


void
BLKSparseMatrixC::reSolve(int nRHS, DComplex *RHS) {

 solveTime -= getTime();

 int i;
 int total = numUncon*nRHS;

 DComplex *s = new DComplex[total];
 DComplex *t = new DComplex[total];

 _FORTRAN(zblkslvp)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp, nRHS, RHS,
                   numUncon, s, t);

 for (i=0; i<total; ++i)
   RHS[i] = s[i];

 delete[] t;
 delete[] s;

 solveTime += getTime();

}


void
BLKSparseMatrixC::reSolve(FullMC *mat) {

 solveTime -= getTime();

 int i, j;
 int nRHS = mat->numCol();
 int length = mat->numRow();

 DComplex *s = new DComplex[length*nRHS];
 DComplex *t = new DComplex[length*nRHS];
 DComplex *a = new DComplex[length*nRHS];

 for (j=0; j<nRHS; ++j)
   for (i=0; i<length; ++i)
     a[i+j*length] = (*mat)[i][j];

 _FORTRAN(zblkslvp)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp, nRHS, a,
                   length, s, t);

 for (j=0; j<nRHS; ++j)
   for (i=0; i<length; ++i)
     (*mat)[i][j] = s[i+j*length];

 delete[] s;
 delete[] t;
 delete[] a;

 solveTime += getTime();

}


void
BLKSparseMatrixC::zeroAll()
{
  int i;
  for(i=0; i < xlnz[numUncon]; ++i)
    lnz[i] = DComplex(0.0, 0.0);
}


void
BLKSparseMatrixC::add(FullSquareMatrix &kel, int *dofs) {

 int i, j, k,rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;

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
             lnz[ position - 2 - offset] += DComplex(kel[i][j], 0.0);
             break;
           }
         }
       }

     }
 }

}


void
BLKSparseMatrixC::add(AssembledFullMC &kel, int *dofs) {

 int i, j, k,rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;

 int kndof = kel.dim();

 for( i = 0; i < kndof; ++i ) {
   if ((rowi = unconstrNum[dofs[i]]) == -1 )
        continue;
   p1     = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper = invsuper[p1] - 1;
   fstcol = xsuper[csuper] - 1;
   lxbeg  = xlindx[csuper]-1;
   lxend  = xlindx[csuper+1]-1;
   for ( j = 0; j < kndof; ++j ) {
       if ((colj = unconstrNum[dofs[j]]) == -1 ) 
          continue;
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


void
BLKSparseMatrixC::addImaginary(FullSquareMatrix &kel, int *dofs)
{

 int i, j, k,rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;

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
             lnz[ position - 2 - offset] += DComplex(0.0, kel[i][j]);
             break;
           }
         }
       }

     }
 }

}


void
BLKSparseMatrixC::add(FullSquareMatrixC &kel, int *dofs)
{

 int i, j, k,rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;

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

// OLD COMMENT FROM REAL BLKSPARSEMATRIX
// This routine has to be converted to input to the LNZ double
// precision array instead of the temporary unonz double precision
// array.

void
BLKSparseMatrixC::add(FullM &knd, int fRow, int fCol)
{
  int offset,k,p1,position, csuper,fstcol, lxbeg,lxend;
  int iCol, iRow;
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

    for(iRow = 0; iRow < nrow; ++iRow) {
       //int rowi = fRow + iRow + 1;
       int rowi = fRow + iRow;
// BEGIN
     p1     = invp[rowi] - 1;
     position = xlnz[p1+1]; 
     csuper = invsuper[p1] - 1;
     fstcol = xsuper[csuper] - 1;
     lxbeg  = xlindx[csuper]-1;
     lxend  = xlindx[csuper+1]-1;
// END
     for(iCol = 0; iCol < ncol; ++iCol) {
          //int colj = fCol + iCol + 1;
          int colj = fCol + iCol ;
// BEGIN
       int irow   = invp[colj] - 1;
       if ( irow >= fstcol ) {
         offset = lxend - lxbeg;
         for(k=lxbeg; k<lxend; ++k) {
           offset -= 1;
           if(lindx[k]-1 == irow) { 
             lnz[ position - 2 - offset] += DComplex(knd[iRow][iCol], 0.0);
             break;
           }
         }
       }
//END
       }
    }
}


void
BLKSparseMatrixC::addImaginary(FullM &knd, int fRow, int fCol)
{
  int offset,k,p1,position, csuper,fstcol, lxbeg,lxend;
  int iCol, iRow;
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

    for(iRow = 0; iRow < nrow; ++iRow) {
       //int rowi = fRow + iRow + 1;
       int rowi = fRow + iRow;
// BEGIN
     p1     = invp[rowi] - 1;
     position = xlnz[p1+1];
     csuper = invsuper[p1] - 1;
     fstcol = xsuper[csuper] - 1;
     lxbeg  = xlindx[csuper]-1;
     lxend  = xlindx[csuper+1]-1;
// END
     for(iCol = 0; iCol < ncol; ++iCol) {
          //int colj = fCol + iCol + 1;
          int colj = fCol + iCol ;
// BEGIN
       int irow   = invp[colj] - 1;
       if ( irow >= fstcol ) {
         offset = lxend - lxbeg;
         for(k=lxbeg; k<lxend; ++k) {
           offset -= 1;
           if(lindx[k]-1 == irow) {
             lnz[ position - 2 - offset] += DComplex(0.0, knd[iRow][iCol]);
             break;
           }
         }
       }
//END
       }
    }
}


void
BLKSparseMatrixC::add(FullMC &knd, int fRow, int fCol)
{
  int offset,k,p1,position, csuper,fstcol, lxbeg,lxend;
  int iCol, iRow;
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

    for(iRow = 0; iRow < nrow; ++iRow) {
       //int rowi = fRow + iRow + 1;
       int rowi = fRow + iRow;
// BEGIN
     p1     = invp[rowi] - 1;
     position = xlnz[p1+1];
     csuper = invsuper[p1] - 1;
     fstcol = xsuper[csuper] - 1;
     lxbeg  = xlindx[csuper]-1;
     lxend  = xlindx[csuper+1]-1;
// END
     for(iCol = 0; iCol < ncol; ++iCol) {
          //int colj = fCol + iCol + 1;
          int colj = fCol + iCol ;
// BEGIN
       int irow   = invp[colj] - 1;
       if ( irow >= fstcol ) {
         offset = lxend - lxbeg;
         for(k=lxbeg; k<lxend; ++k) {
           offset -= 1;
           if(lindx[k]-1 == irow) {
             lnz[ position - 2 - offset] += knd[iRow][iCol];
             break;
           }
         }
       }
//END
       }
    }
}


void
BLKSparseMatrixC::printAll()
{
  fprintf(stderr, " numUncon       = %d\n",numUncon);
  fprintf(stderr, " neq            = %d\n",neq);
  fprintf(stderr, " numConstrained = %d\n",numConstrained);
}


DComplex
BLKSparseMatrixC::diagComplex(int dof)
{
 int k,colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 p1 = invp[dof]-1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 colj = dof;
 int irow   = invp[colj] - 1;
 if ( irow >= fstcol ) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
      offset -= 1;
      if(lindx[k]-1 == irow) 
            return lnz[ position - 2 - offset];
   }
 }
 return DComplex(0.0, 0.0);
}


BLKSparseMatrixC::BLKSparseMatrixC(Connectivity *cn, DofSetArray *_dsa, 
                                 int *glInternalMap, double _tol) :
  SparseData(_dsa, glInternalMap, cn, 1)
{
  numrbm = 0;
  tol    = _tol;
  allocateMemory();
}


BLKSparseMatrixC::BLKSparseMatrixC(Connectivity *cn,EqNumberer *_dsa,
                  double _tol) : SparseData(cn,_dsa,_tol) {

  numrbm = 0;
  tol    = _tol;
  allocateMemory();

}


void
BLKSparseMatrixC::allocateMemory() {

  solveTime = 0.0;

  // Rigid body mode information
  numrbm = 0;

  snode  = new int[numUncon];
  xsuper = new int[numUncon+1];
  xlindx = new int[numUncon+1];
  xlnz   = new int[numUncon+1];
  iwork  = new int[7*numUncon+3];
  perm   = new int[numUncon];
  invp   = new int[numUncon];
  colcnt = new int[numUncon];

  lindx  = new int[xunonz[numUncon]];

  int nnza = xunonz[numUncon] - 1 - numUncon;

  int i,j,k,iflag;

  // maxsup = maximum number super nodes (parameter)

  int maxsup = 150;
  //int maxsup = 125;
  //int maxsup = 100;
  //int maxsup =  75;
  //int maxsup =  50;
  //int maxsup =  40;
  //int maxsup =  35;
  //int maxsup =  30;
  //int maxsup =  25;
  //int maxsup =  20;
  //int maxsup =  10;

  adj = new int[nnza];

  for (k=0, j=0; k< numUncon; k++)
     for (i=xunonz[k]; i <=xunonz[k+1]-1; i++)
        if ( rowu[i-1] == k+1 );
        else  {
          adj[j++] = rowu[i-1];
        }

  xadj = new int[numUncon + 1];

  for (i=0; i< numUncon+1; i++)
     xadj[i] = xunonz[i] - i ;

  //DEBUG
  //fprintf(stderr,"xunonz\n");
  //for (i=0; i< numUncon+1; i++)
  //  fprintf(stderr,"%d\n",xunonz[i]);
  //DEBUG
  //fprintf(stderr,"rowu\n");
  //for (k=0; k < numUncon; k++)
  //     for (int j=xunonz[k]-1; j < xunonz[k+1]-1; j++)
  //         fprintf(stderr,"%d\n",rowu[j]);
  
  delete [] xunonz; xunonz = 0;
  delete [] rowu;   rowu   = 0;

//-------------------------------------------
// Copy matrix structure from (XADJ,ADJ) to
// (XLINDX,LINDX) (because matrix structure is
// destroyed by the minimum degree ordering
// subroutine).
// -------------------------------------------

  for (k=0; k< numUncon+1; k++)
     xlindx[k] = xadj[k];

  for (k=0; k<nnza; k++)
     lindx[k] = adj[k];

// -------------------------------------------------------
// STEP 1:   ORDMMD ...  multiple minimum degree ordering.
//
//           Input:      NUMUNCON, XLINDX, LINDX, IWSIZE
//           Output:     PERM, INVP, NSUB, IFLAG
//           Work:       IWORK(4*N)
// -------------------------------------------------------

  iwsiz = 4 * numUncon;

  _FORTRAN(ordmmd2)(numUncon,  xlindx,   lindx, invp, perm,
                      iwsiz, iwork,  nsub, iflag);

//               IFLAG =  0: successful ordering.
//               IFLAG = 11: insufficient work space in IWORK.

  if(iflag == 11)
    fprintf(stderr,"*** ERROR: not enough space for iwork temperary array\n");

  iwsiz = 7 * numUncon + 3;

  // Rank deficiency information
  // defblk = size of last block to perform full pivoting on
  defblk = 0;

//      ***********************
//      Symbolic factorization.
//      ***********************
//      -------------------------------------------------------------
//      STEP 2:
//      SFINIT ...  symbolic factorization initialization, which
//                   computes supernode partition and storage
//                  requirements for symbolic factorization;
//                  new ordering is a postordering of the nodal
//                   elimination tree.
//
//       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, MAXSUP, DEFBLK,
//                  IWSIZE
//       Output:     PERM, INVP, COLCNT, NNZL, NSUB, NSUPER, XSUPER,
//                   SNODE, IFLAG
//       Work:       IWORK(7*N+3) ... the max required any subroutine.
//       -------------------------------------------------------------


  _FORTRAN(sfinit)(numUncon, nnza,   xadj,   adj,   perm,
                   invp,     maxsup, defblk, colcnt, nnzl,
                   nsub,     nsuper,  xsuper, snode, iwsiz,
                   iwork,    iflag);
//IFLAG =  0: successful symbolic factorization
//            initialization.
//IFLAG = 21: insufficent work space in IWORK.
  if (iflag == 21) {
     fprintf(stderr," *** ERROR: Insufficient working storage for SFINIT\n");
     exit(1);
  }

//       ------------------------------------------------------
//       STEP 3:
//       SYMFCT ...  perform supernodal symbolic factorization.
//
//       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, COLCNT,
//                   NSUPER, XSUPER, SNODE , NSUB, IWSIZE
//       Output:     XLINDX, LINDX, XLNZ, IFLAG
//       WORK:       IWORK(NSUPER+2*N+1)
//
//       No longer needed: ADJ, XADJ, COLCNT
//       ------------------------------------------------------

  iwsiz = nsuper + 2 * numUncon + 1;

  _FORTRAN(symfct)(numUncon, nnza,   xadj,   adj,    perm,
                   invp,     colcnt, nsuper, xsuper, snode,
                   nsub,     xlindx, lindx, xlnz,    iwsiz,
                   iwork,    iflag);

//               IFLAG =  0: no error.
//               IFLAG = 22: insufficient work space in IWORK.
//               IFLAG = 23: inconsistancy in the input.

  // ERROR Checking
  if (iflag == 22) {
     fprintf(stderr," *** ERROR: (SYMFCT) Memory for iwork is too small\n");
     exit(1);
  } else if (iflag == 23) {
     fprintf(stderr," *** ERROR: (SYMFCT) Inconsistency in input\n");
     exit(1);
  }

  // Delete memory that is not used any longer
  delete [] xadj;   xadj   = 0;
  delete [] adj;    adj    = 0;
  delete [] colcnt; colcnt = 0;

  dbg_alloca(0);

  lnz = new DComplex[xlnz[numUncon]];

  // zero lnz and create invsuper array
  invsuper = new int[numUncon];

  _FORTRAN(zzeromat)(numUncon,nsuper,xsuper,xlnz,lnz,invsuper);

}

