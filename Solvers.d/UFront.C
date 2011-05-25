#include <cstdio>
#include <cstdlib>

#ifdef sgi
#include <osfcn.h>
#include <ieeefp.h>
#endif

#include <Solvers.d/UFront.h>
#include <Utils.d/dofset.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/mathUtility.h>
#include <Timers.d/GetTime.h>


#define FREE  -1
#define FIXED -2
#define NEW   -3

// nops = number of operations
double nops = 0;


void
rankNUpdate_3(double *k1, double *k2, double *k3, int len, double a[6][3], double (*r)[6])
{
     double a11, a12, a13, a14, a15, a16;
     double a21, a22, a23, a24, a25, a26;
     double a31, a32, a33, a34, a35, a36;

     a11 = a[0][0]; a12 = a[1][0]; a13 = a[2][0]; a14 = a[3][0]; a15 = a[4][0]; a16 = a[5][0];
     a21 = a[0][1]; a22 = a[1][1]; a23 = a[2][1]; a24 = a[3][1]; a25 = a[4][1]; a26 = a[5][1];
     a31 = a[0][2]; a32 = a[1][2]; a33 = a[2][2]; a34 = a[3][2]; a35 = a[4][2]; a36 = a[5][2];

 int i;
 #pragma ivdep
 for(i= 0; i < len; ++i) {
      k1[i] += a11*r[i][0] + a12*r[i][1] + a13*r[i][2] +
               a14*r[i][3] + a15*r[i][4] + a16*r[i][5];
      k2[i] += a21*r[i][0] + a22*r[i][1] + a23*r[i][2] +
               a24*r[i][3] + a25*r[i][4] + a26*r[i][5];
      k3[i] += a31*r[i][0] + a32*r[i][1] + a33*r[i][2] +
               a34*r[i][3] + a35*r[i][4] + a36*r[i][5];
 }
}

void
rankNUpdate_2(double *k1, double *k2, int len, double a[6][3], double (*r)[6])
{
     double a11, a12, a13, a14, a15, a16;
     double a21, a22, a23, a24, a25, a26;

     a11 = a[0][0]; a12 = a[1][0]; a13 = a[2][0]; a14 = a[3][0]; a15 = a[4][0]; a16 = a[5][0];
     a21 = a[0][1]; a22 = a[1][1]; a23 = a[2][1]; a24 = a[3][1]; a25 = a[4][1]; a26 = a[5][1];

 int i;
 #pragma ivdep
 for(i= 0; i < len; ++i) {
      k1[i] += a11*r[i][0] + a12*r[i][1] + a13*r[i][2] +
               a14*r[i][3] + a15*r[i][4] + a16*r[i][5];
      k2[i] += a21*r[i][0] + a22*r[i][1] + a23*r[i][2] +
               a24*r[i][3] + a25*r[i][4] + a26*r[i][5];
 }
}



void zeroinit(double *d, int i)
{
 while(i--)
  d[i] = 0.0;
}

UFront::UFront(int _ndof, int _maxfrsize, double prec, Rbm *_rbm, int)
 {
  solveTime = 0.0;
  rbm = _rbm;
  if(rbm)
    numrbm = rbm->numRBM();
  else
    numrbm = 0;

  // Compute memory used by Frontal Solver
  memUsed = -memoryUsed();

  int i, j, id, iloc;
  ndof = _ndof;
  precision = prec;
  maxfrsize = _maxfrsize;
  loc = new int[ndof];
  for(iloc = 0; iloc < ndof ; ++iloc)
    loc[iloc] = NEW;
  invloc = new int[maxfrsize] ;
  for(i = 0; i < maxfrsize; ++i)
    invloc[i] = FREE ;

  data = new double[(maxfrsize*(maxfrsize+1))/2 + 2*maxfrsize] ;
  if(data == 0) exit(-1) ;   // Failed to allocate data space

  k = new double*[maxfrsize] ;
  /* Initialize k */
  for(j = 0, id=0; j < maxfrsize; ++j) {
    k[j] = data+id;
    id += j+3;
  }

  rhs = new double[maxfrsize] ;

  zeroinit(data,(maxfrsize*(maxfrsize+1))/2);
  zeroinit(rhs,maxfrsize);

  int maxSaved = (ndof/UNROLL) +1;
  currentFrontSize = 0;
  nSaved = 0;
  row = new double[maxfrsize][UNROLL];
  column = new double[maxfrsize][UNROLL];
  order = new int[maxSaved][UNROLL];
  rowlen = new int[maxSaved];
  savedFrontPos = new int[maxSaved][UNROLL];
  savedDiagBlock = new double[maxSaved][UNROLL][UNROLL];
  savedRHS = new double[maxSaved][UNROLL];
  savedRows = new double *[maxSaved][UNROLL];
  nRows = 0;
  dsa = 0;
}

UFront::UFront(ConstrainedDSA* c_dsa, int _maxfrsize, double prec, Rbm *_rbm, int)
{
  solveTime = 0.0;
  // Compute memory used by Frontal Solver
  memUsed = -memoryUsed();

  rbm = _rbm;
  if(rbm)
    numrbm = rbm->numRBM();
  else
    numrbm = 0;

  int i, j, id, iloc;
  ndof = c_dsa->size();
  precision = prec;
  maxfrsize = _maxfrsize;
  loc = new int[ndof];
  for(iloc = 0; iloc < ndof ; ++iloc)
    loc[iloc] = NEW;
  invloc = new int[maxfrsize];
  for(i = 0; i < maxfrsize; ++i)
    invloc[i] = FREE ;

  data = new double[(maxfrsize*(maxfrsize+1))/2 + 2*maxfrsize];
  if(data == 0) exit(-1);   // Failed to allocate data space

  k = new double*[maxfrsize] ;
  /* Initialize k */
  for(j = 0, id=0; j < maxfrsize; ++j) {
    k[j] = data+id ;
    id += j+3;
  }

  rhs = new double[maxfrsize];

  zeroinit(data,(maxfrsize*(maxfrsize+1))/2);
  zeroinit(rhs,maxfrsize);

  int maxSaved = (ndof/UNROLL) +1;
  currentFrontSize = 0;
  nSaved = 0;
  row = new double[maxfrsize][UNROLL];
  column = new double[maxfrsize][UNROLL];
  order = new int[maxSaved][UNROLL];
  rowlen = new int[maxSaved];
  savedFrontPos = new int[maxSaved][UNROLL];
  savedDiagBlock = new double[maxSaved][UNROLL][UNROLL];
  savedRHS = new double[maxSaved][UNROLL];
  savedRows = new double *[maxSaved][UNROLL];
  nRows = 0;
  dsa = c_dsa;
}

UFront::~UFront()
{
  if(loc) { delete [] loc; loc = 0; }
  if(invloc) { delete [] invloc; invloc = 0; }
  if(data) { delete [] data; data = 0; }
  if(k) { delete [] k; k = 0; }
  if(rhs) { delete [] rhs; rhs = 0; }
  if(row) { delete [] row; row = 0; }
  if(column) { delete [] column; column = 0; }
  if(order) { delete [] order; order = 0; }
  if(rowlen) { delete [] rowlen; rowlen = 0; }
  if(savedFrontPos) { delete [] savedFrontPos; savedFrontPos = 0; }
  if(savedDiagBlock) { delete [] savedDiagBlock; savedDiagBlock = 0; }
  if(savedRHS) { delete [] savedRHS; savedRHS = 0; }
  if(savedRows) { delete [] savedRows; savedRows = 0; }
}

int UFront::locate(int i)
 {
  int j,l;
  if(loc[i] != NEW) return loc[i] ;
  for(j = 0 ; j < maxfrsize ; ++ j)
    if(invloc[j] < 0)
      {
        invloc[j] = i;
        loc[i] = j ;
        if (j == currentFrontSize) {
           for(l = 0; l < UNROLL; ++l)
             column[j][l] = row[j][l] = 0.0;
           ++currentFrontSize ;
        }
        return j ;
      }
  fprintf(stderr," Front was too small!!!\n") ;
  exit(-1) ;
  return -1 ;
 }

void UFront::mark_fixed(int dof)
{
 loc[dof] = FIXED ;
}


void UFront::addkel(SymFullMatrix &Kel, int *locel, double *q)
{
  int i,j,loci,locj;
  for(i = 0; i < Kel.dim() ; ++i) {
    int ii = locel[i];  //  ii is the global number of the DOF
    if(ii < 0) continue;
    loci = locate(ii);  //  loci is the local (to the front) number
    for(j = 0; j <= i; ++j) {
      int jj = locel[j];
      if(jj < 0) continue;
      locj = locate(jj);
      if(loci == FIXED) {
        if(locj == FIXED)
          continue;   // Both are fixed, => skip
        else
          rhs[locj] -= q[ii]*Kel[i][j];   // Only ii is fixed
      }
      else {  // ii is NOT Fixed
        if(locj == FIXED)
          rhs[loci] -= q[jj]*Kel[i][j];
        else {
          if(loci <= locj)
            k[locj][loci] += Kel[i][j];
          else
            k[loci][locj] += Kel[i][j];
          if(ii == jj && i != j)
            k[loci][locj] += Kel[i][j];
         }
       } // if(ii == FIXED)
     }
   }
}

void UFront::addkel(FullSquareMatrix &Kel, int *locel, double *q)
{
  int i,j,loci,locj;

  for(i = 0; i < Kel.dim() ; ++i) {
    int ii = locel[i];  //  ii is the global number of the DOF
    if(dsa) ii = dsa->getRCN(ii); // Constrained numbering
    if(ii < 0) continue;
    loci = locate(ii);  //  loci is the local (to the front) number

    for(j = 0; j <= i; ++j) {
      int jj = locel[j];
      if(dsa) jj = dsa->getRCN(jj); // Constrained numbering
      if(jj < 0) continue;
      locj = locate(jj);
      if(loci == FIXED) {
        if(locj == FIXED)
          continue;   // Both are fixed, => skip
        else
          if(q) rhs[locj] -= q[ii]*Kel[i][j];   // Only ii is fixed
      }
      else { // ii is NOT Fixed
        if(locj == FIXED)
          rhs[loci] -= q[jj]*Kel[i][j];
        else {
          if(loci <= locj)
            k[locj][loci] += Kel[i][j];
          else
            k[loci][locj] += Kel[i][j];
          if(ii == jj && i != j)
            k[loci][locj] += Kel[i][j];
         }
       } // if(ii == FIXED)
     }
   }
}


void UFront::addload(int dof, double v)
{
 int locdof = locate(dof);
 if(locdof >= 0)
   rhs[locdof] += v;
}

void UFront::elim(int dof)
{
 if(loc[dof] == FIXED) return;
 collect(dof);
 
 if(nRows == UNROLL) {
   subFactor();
   rankNUpdate();
   saveEliminatedRows();
 }
}

void UFront::addToDiag(int dof, double v)
{
 int loci = loc[dof];
 if(loci == FIXED) return;
 k[loci][loci] += v;
}

void UFront::collect(int dof)
{
 int locdof = locate(dof);
 
 int i, j;
 for(j = 0; j < locdof; ++j)
  {
   row[j][nRows] = k[locdof][j];
   k[locdof][j] = 0.0;
  } 
 diagBlock[nRows][nRows] = k[locdof][locdof];
 k[locdof][locdof] = 0.0;
// Due to unrolling:
 k[locdof][locdof+1] = 0.0;
 k[locdof][locdof+2] = 0.0;

 for(i = locdof+1; i < currentFrontSize; ++i)
  {
   row[i][nRows] = k[i][locdof];
   k[i][locdof] = 0.0;
  }

 row[locdof][nRows] = 0.0;
 factRHS[nRows] = rhs[locdof];
 rhs[locdof] = 0.0;
 rowFrontPos[nRows] = locdof;

 for(i = 0; i < nRows; ++i) {
   diagBlock[nRows][i] = diagBlock[i][nRows] = row[locdof][i];
   row[locdof][i] = 0.0;
 }

 invloc[locdof] = FREE;
 if(currentFrontSize == locdof)
  {
   while(currentFrontSize > 0 && invloc[currentFrontSize] == FREE)
     currentFrontSize--;
  }
 rowFrontPos[nRows] = locdof;
 rowDOF[nRows] = dof;
 nRows++;
}

void UFront::subFactor()
{
 int iRow, jRow, j;
 // Save the rows as columns before subelimination
 for(iRow = 0; iRow < nRows; ++ iRow)
   for(j = 0; j < currentFrontSize; ++j)
     column[j][iRow] = row[j][iRow];
 // Perform pivoting
 for(iRow = 0; iRow < nRows; ++ iRow) {
    // Find the maximum diagBlock[iRow][jRow]
    int maxDIndex = iRow;
    double maxD = diagBlock[iRow][iRow];
    /*for(jRow = iRow+1; jRow < nRows; ++jRow)
       if(abs(diagBlock[jRow][jRow]) > maxD) {
          maxD = abs(diagBlock[jRow][jRow]);
          maxDIndex = jRow;
        }
    */
    // Test singularity
    if(maxD < precision)
    {
     fprintf(stderr,"DOF is singular %e %e\n",precision, maxD);
     diagBlock[iRow][iRow] = 0.0;
     continue;
    }
    // Swap all the data
    if(maxDIndex != iRow) { // Swap lines
      for(j=0; j < nRows; ++j) {
        double tmp = diagBlock[iRow][j];
        diagBlock[iRow][j] = diagBlock[maxDIndex][j];
        diagBlock[maxDIndex][j] = tmp;
      }
      for(j=0; j < nRows; ++j) {
        double tmp = diagBlock[j][iRow];
        diagBlock[j][iRow] = diagBlock[j][maxDIndex];
        diagBlock[j][maxDIndex] = tmp;
      }
      // fprintf(stderr,"Pivot %d %d %e %e\n",maxDIndex , iRow, diagBlock[iRow][iRow], maxD);
      int tmp = rowDOF[iRow];
      rowDOF[iRow] = rowDOF[maxDIndex];
      rowDOF[maxDIndex] = tmp;
      tmp = rowFrontPos[iRow];
      rowFrontPos[iRow] = rowFrontPos[maxDIndex];
      rowFrontPos[maxDIndex] = tmp;

/*
      double *tCol = column[iRow];
      column[iRow] = column[maxDIndex];
      column[maxDIndex] = tCol;
      tCol = row[iRow];
      row[iRow] = row[maxDIndex];
      row[maxDIndex] = tCol;
*/
      for(j = 0; j < currentFrontSize; ++j) {
        double tmp = column[j][iRow];
        column[j][iRow] = column[j][maxDIndex];
        column[j][maxDIndex]=tmp;
        tmp = row[j][iRow];
        row[j][iRow] = row[j][maxDIndex];
        row[j][maxDIndex]=tmp;
      }
    }

  // subelimination
  for(jRow = iRow+1; jRow < nRows; ++jRow) {
    double coef = diagBlock[iRow][jRow] / diagBlock[iRow][iRow];
    for(j = iRow; j < nRows; ++j)
      diagBlock[jRow][j] -= coef * diagBlock[iRow][j];
    for(j = 0; j < currentFrontSize; ++j)
      row[j][jRow] -= coef*row[j][iRow];
    factRHS[jRow] -= coef*factRHS[iRow];
  }
 }
}

void
UFront::rankNUpdate() {
  double *ir[3];
  int i,j,l,m;
  for(i = 0; i < currentFrontSize; i += 3) {
     ir[0] = ir[1] = ir[2] = 0;
     double a[6][3];
     double a11, a12, a13, a14, a15, a16;
     int nv =0;
     for(j =0; j < 3; ++j) {
       if(i+j >= currentFrontSize || invloc[i+j] == FREE)
          continue;
       ir[nv] = k[i+j];
       for(l=0; l < UNROLL; ++l)
         a[l][nv] = column[i+j][l];
       for(l = 0; l < UNROLL; ++l) {
         double coef = (diagBlock[l][l] != 0.0) ? 
                        a[l][nv] = -a[l][nv]/diagBlock[l][l] : 0;
         for(m = l+1; m < UNROLL; ++m)
           a[m][nv] += coef*diagBlock[l][m];
       }
       for(l = 0; l < UNROLL; ++l) 
          rhs[i+j] += a[l][nv]*factRHS[l];
       nv++;
     }
     if(nv == 0) continue;
     if(nv == 1) {
     a11 = a[0][0]; a12 = a[1][0]; a13 = a[2][0]; a14 = a[3][0]; a15 = a[4][0]; a16 = a[5][0];
       for(j = 0; j < i+3; ++j) {
         ir[0][j] += a11*row[j][0] + a12*row[j][1] + a13*row[j][2]
                   + a14*row[j][3] + a15*row[j][4] + a16*row[j][5];
       }
       nops += 12*(i+3);
     }
     if(nv == 2) {
     /*a11 = a[0][0]; a12 = a[1][0]; a13 = a[2][0]; a14 = a[3][0]; a15 = a[4][0]; a16 = a[5][0];
     a21 = a[0][1]; a22 = a[1][1]; a23 = a[2][1]; a24 = a[3][1]; a25 = a[4][1]; a26 = a[5][1];
       for(j = 0; j < i+3; ++j) {
         ir[0][j] += a11*row[j][0] + a12*row[j][1] + a13*row[j][2]
                   + a14*row[j][3] + a15*row[j][4] + a16*row[j][5];
         ir[1][j] += a21*row[j][0] + a22*row[j][1] + a23*row[j][2]
                   + a24*row[j][3] + a25*row[j][4] + a26*row[j][5];
       }*/
       ::rankNUpdate_2(ir[0],ir[1], i+3, a, row);
       nops += 24*(i+3);
     }
     if(nv == 3) {
/*
     a11 = a[0][0]; a12 = a[1][0]; a13 = a[2][0]; a14 = a[3][0]; a15 = a[4][0]; a16 = a[5][0];
     a21 = a[0][1]; a22 = a[1][1]; a23 = a[2][1]; a24 = a[3][1]; a25 = a[4][1]; a26 = a[5][1];
     a31 = a[0][2]; a32 = a[1][2]; a33 = a[2][2]; a34 = a[3][2]; a35 = a[4][2]; a36 = a[5][2];

#pragma ivdep
       for(j = 0; j < i+3; ++j) {
      ir[0][j] += a11*row[j][0] + a12*row[j][1] + a13*row[j][2] +
                  a14*row[j][3] + a15*row[j][4] + a16*row[j][5];
      ir[1][j] += a21*row[j][0] + a22*row[j][1] + a23*row[j][2] +
                  a24*row[j][3] + a25*row[j][4] + a26*row[j][5];
      ir[2][j] += a31*row[j][0] + a32*row[j][1] + a33*row[j][2] +
                  a34*row[j][3] + a35*row[j][4] + a36*row[j][5];
       }*/
       ::rankNUpdate_3(ir[0],ir[1],ir[2], i+3, a, row);
       nops += 36*(i+3);
     }
  }
}


void UFront::saveEliminatedRows()
{
 int i,j;
 for(i = 0; i < nRows; ++i) {
   order[nSaved][i] = rowDOF[i];
   for(j= 0; j < nRows; ++j)
      savedDiagBlock[nSaved][i][j] = diagBlock[i][j];
   
   savedRHS[nSaved][i] = factRHS[i];
   if(currentFrontSize > 0)
     savedRows[nSaved][i] = new double[currentFrontSize];
   for(j =0; j < currentFrontSize ; ++j)
     savedRows[nSaved][i][j] = row[j][i];
   rowlen[nSaved] = currentFrontSize;
   savedFrontPos[nSaved][i] = rowFrontPos[i];
 }
 for(; i < UNROLL; ++i) {
   savedFrontPos[nSaved][i] = -1;
   order[nSaved][i] = -1;
 }
 nSaved++;
 nRows = 0;
}

void UFront::backsub(double *q)
{
 double *qfr = new double[maxfrsize];
 int nsv = nSaved-1;
 int i,j;
 // First do a backsub just on the last one
 for(i = nRows; i-- > 0; ) {
    for(j = i+1; j < nRows; ++j) {
      savedRHS[nsv][i] -= savedDiagBlock[nsv][i][j] * savedRHS[nsv][j];
    }
    if(abs(savedDiagBlock[nsv][i][i]) < precision)
       savedRHS[nsv][i]  = 0;
     else
       savedRHS[nsv][i] /= savedDiagBlock[nsv][i][i];
    qfr[savedFrontPos[nsv][i]] = q[order[nsv][i]] = savedRHS[nsv][i];
 }
 while(nsv-- > 0) {
   for(i = 0; i < rowlen[nsv]; ++ i) {
     for(j = 0; j < UNROLL; ++j)
       savedRHS[nsv][j] -= savedRows[nsv][j][i]*qfr[i];
   }
   for(i = UNROLL; i-- > 0; ) {
     if(abs(savedDiagBlock[nsv][i][i]) < precision)
       {
        savedRHS[nsv][i] = qfr[savedFrontPos[nsv][i]] = q[order[nsv][i]] = 0.0;
        continue;
       }
     for(j = i+1; j < UNROLL; ++j) 
        savedRHS[nsv][i] -= savedDiagBlock[nsv][i][j] * savedRHS[nsv][j];

     savedRHS[nsv][i] /= savedDiagBlock[nsv][i][i];
   }
   for(i = UNROLL; i-- > 0; )
     q[order[nsv][i]] = qfr[savedFrontPos[nsv][i]] = savedRHS[nsv][i]; 
 }
}

void
UFront::finishUpdate() {
 subFactor();
 int t = nRows;
 saveEliminatedRows();
 nRows = t; // Not clean should modify backsub and not call save Eliminated
}

void
UFront::solve(Vector &rightHandSide, Vector &solution)
{
 fprintf(stderr," ... Made %f Million Floating Points\n",nops/1e6);
 solution = rightHandSide;
 reSolve(solution.data());
}

void
UFront::solve(double *rightHandSide, double* solution)
{
 int i;
 for(i=0; i<neqs(); ++i)
   solution[i] = rightHandSide[i];

 reSolve(solution);
}

void
UFront::reSolve(double *rightHandSide)
{
 solveTime -= getTime();

 double *qfr = new double[maxfrsize];
 zeroinit(qfr,maxfrsize);
 int i,j,k;

 // Forward substitute
 for(i = 0; i < nSaved; ++i) {
   int nUnroll = (i < nSaved-1) ?  UNROLL : nRows;
   for(j =0; j < nUnroll; ++j)
     factRHS[j] = 0.0;
   for(j =0; j < nUnroll; ++j) {
     if(abs(savedDiagBlock[i][j][j])  < precision)
      {
        // Should do something
        continue;
      }
     double coef = 1.0/savedDiagBlock[i][j][j];
     factRHS[j] += rightHandSide[order[i][j]];
     factRHS[j] += qfr[savedFrontPos[i][j]];
     qfr[savedFrontPos[i][j]] = 0.0;
     for(k = j+1; k < nUnroll; ++k)
        factRHS[k] -= coef*savedDiagBlock[i][j][k]*factRHS[j];
      for(k = 0; k < rowlen[i] ; ++k)
       {
          qfr[k] -= coef*savedRows[i][j][k]*factRHS[j];
       }
    }
   for(j =0; j < nUnroll; ++j)
     rightHandSide[order[i][j]] = factRHS[j];
  }
 int nsv = nSaved-1;
 // First do a backsub just on the last one
 for(j =0; j < nRows; ++j)
   factRHS[j] = rightHandSide[order[nsv][j]];
 for(i = nRows; i-- > 0; ) {
    for(j = i+1; j < nRows; ++j) {
      factRHS[i] -= savedDiagBlock[nsv][i][j] * factRHS[j];
    }
    if(abs(savedDiagBlock[nsv][i][i]) < precision)
        rightHandSide[order[nsv][i]] = factRHS[i]  = 0;
     else
        rightHandSide[order[nsv][i]] = (factRHS[i] /= savedDiagBlock[nsv][i][i]);
    qfr[savedFrontPos[nsv][i]] = rightHandSide[order[nsv][i]];
 }
 while(nsv-- > 0) {
   for(j = 0; j < UNROLL; ++j)
      factRHS[j] = rightHandSide[order[nsv][j]];

   for(i = 0; i < rowlen[nsv]; ++ i) {
     for(j = 0; j < UNROLL; ++j)
       factRHS[j] -= savedRows[nsv][j][i]*qfr[i];
   }
   for(i = UNROLL; i-- > 0; ) {
     if(abs(savedDiagBlock[nsv][i][i]) < precision)
       {
        factRHS[i] = 
          qfr[savedFrontPos[nsv][i]] = rightHandSide[order[nsv][i]] = 0.0;
        continue;
       }
     for(j = i+1; j < UNROLL; ++j)
        factRHS[i] -= savedDiagBlock[nsv][i][j] * factRHS[j];

     factRHS[i] /= savedDiagBlock[nsv][i][i];
   }
   for(i = UNROLL; i-- > 0; )
     rightHandSide[order[nsv][i]] = qfr[savedFrontPos[nsv][i]] = factRHS[i];
 }

 solveTime += getTime();
 memUsed += memoryUsed();
}


void
UFront::getRBMs(double *rigidBodyModes)
{ 
  rbm->getRBMs(rigidBodyModes);
}

void
UFront::getRBMs(Vector *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

void
UFront::getRBMs(VectorSet& rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}
