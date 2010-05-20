#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/FullSquareMatrix.h>
#include <Comm.d/Communicator.h>
#include <Threads.d/Paral.h>

extern double getTime();

#include<Driver.d/Domain.h>
extern Domain *domain;

#define DEBUG_SPOOLES 0  // 1 = print debug stats for factorization, 2 = print for solve also


/*************************************************************************************/
/*************************************************************************************/
#ifdef USE_SPOOLES

extern "C" {
 #include <MT/spoolesMT.h>
}

inline double DVTrace( DV &d )
{
 double t = 0.0;
 for(int i = 0; i < d.size; ++i)
   t += ((i%11)+1)*d.vec[i];
 return t;
}

inline void DenseMtx_setRealOrComplexEntry(DenseMtx *mtx, int row, int col, double value)
{ DenseMtx_setRealEntry(mtx, row, col, value); }

inline void DenseMtx_setRealOrComplexEntry(DenseMtx *mtx, int row, int col, DComplex value)
{ DenseMtx_setComplexEntry(mtx, row, col, value.real(), value.imag()); }

inline void InpMtx_inputRealOrComplexEntry(InpMtx *inpmtx, int row, int col, double value)
{ InpMtx_inputRealEntry(inpmtx, row, col, value); }

inline void InpMtx_inputRealOrComplexEntry(InpMtx *inpmtx, int row, int col, DComplex value)
{ InpMtx_inputComplexEntry(inpmtx, row, col, value.real(), value.imag()); }

inline void DenseMtx_realOrComplexEntries(DenseMtx *mtx, double *entries, int neq)
{ 
  double *temp = DenseMtx_entries(mtx);
  for(int i=0; i<neq; ++i) entries[i] = temp[i];
}


inline void DenseMtx_realOrComplexEntries(DenseMtx *mtx, DComplex *entries, int neq)
{
// RT: fixed an outrageous bug here
  int i;
  for(i=0; i<neq; ++i) {
    int irow = i;
    int jcol = 0;
    double re, im;
    DenseMtx_complexEntry(mtx, irow, jcol, &re, &im);
    DComplex p(re, im);
    entries[i] = p;
  }
}   

template<class Scalar>
class SpoolesType {
};

template <>
class SpoolesType<double> {
 public:
  static const int type = SPOOLES_REAL;
};

template <>
class SpoolesType<complex<double> > {
 public:
  static const int type = SPOOLES_COMPLEX;
};

#endif

template<class Scalar>
GenSpoolesSolver<Scalar>::GenSpoolesSolver(Connectivity *nToN, EqNumberer *_dsa, int *map)
 : SparseData(_dsa, nToN, map)
{
  // constructor for feti-dp coarse problem
  init();
  neq = numUncon; 
  //int nNodes = nToN->csize();
  nNonZero = xunonz[numUncon]-1; 

  unonz = new Scalar[xunonz[numUncon]];
  for(int i = 0; i < nNonZero; ++i)
     unonz[i] = 0.0;

#ifdef USE_SPOOLES
  inpMtx = InpMtx_new();

  // We should know what the number of entries is...
  InpMtx_init(inpMtx, INPMTX_BY_COLUMNS, SpoolesType<Scalar>::type, nNonZero, neq);
#endif
}

template<class Scalar>
GenSpoolesSolver<Scalar>::GenSpoolesSolver(Connectivity *nToN, DofSetArray *_dsa,
                                           ConstrainedDSA *c_dsa, Rbm *_rbm) 
 : SparseData(_dsa,c_dsa,nToN)
{
  init();
  neq = c_dsa->size();
  myMem = 0; // this means don't delete unconstrNum
  unconstrNum = c_dsa->getUnconstrNum();

  //int nNodes = nToN->csize();
  nNonZero = xunonz[numUncon]-1; 

  unonz = new Scalar[xunonz[numUncon]];
  for(int i = 0; i < nNonZero; ++i)
     unonz[i] = 0.0;

#ifdef USE_SPOOLES
  inpMtx = InpMtx_new();

  InpMtx_init(inpMtx, INPMTX_BY_COLUMNS, SpoolesType<Scalar>::type, nNonZero, neq);
#endif

  // PJSA 4-13-05: store pointer to geometric RBM
  if(_rbm) {
    rbm = _rbm;
    nrbm  = rbm->numRBM();
  }
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::add(FullSquareMatrix &kel, int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue; 
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        // if(rowu[m-1] > unconstrNum[dofs[j]]+1) 
        //   fprintf(stderr, "Bigger: %d %d\n", rowu[m-1]-1, unconstrNum[dofs[j]]);
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

/* NOT TESTED
template<class Scalar>
void
GenSpoolesSolver<Scalar>::addBoeing(int nlines, const int *Kai, const int *Kaj,
                                    const double *nz, int *map, Scalar multiplier)
{
  int i, j, m, mstart, mstop;

  for(i = 0; i < nlines; ++i) {
    if(map[i] == -1) continue;
    if(unconstrNum[map[i]] == -1) continue; 
    for(j = Kai[i]; j < Kai[i+1]; ++j) {
      if(unconstrNum[map[Kaj[j-1]-1]] == -1) continue;
      if(unconstrNum[map[Kaj[j-1]-1]] < unconstrNum[map[i]]) continue;
      mstart = xunonz[unconstrNum[map[Kaj[j-1]-1]]];
      mstop  = xunonz[unconstrNum[map[Kaj[j-1]-1]]+1];
      for(m=mstart; m<mstop; ++m) {
        if(rowu[m-1] == (unconstrNum[map[i]] + 1)) {
          unonz[m-1] += (nz[j-1]*multiplier);
          break;
        }
      }
    }
  }
}
*/


template<class Scalar>
void
GenSpoolesSolver<Scalar>::add(int dofi, int dofj, Scalar d)
{
  // WARNING: this adds only [i,j] with  j >= i (i.e. upper part)
  if((dofi < 0) || (dofj < 0)) return;
  int m, mstart, mstop, rowi, colj;
  if(unconstrNum) {
    if((rowi = unconstrNum[dofi]) == -1 || (colj = unconstrNum[dofj]) == -1) return;
  }
  else { rowi = dofi; colj = dofj; }

  if(colj<rowi) { // swap row & col to be in the upper part
    int tmp = colj;
    colj = rowi;
    rowi = tmp;
  }
  // upper part
  mstart = xunonz[colj];
  mstop  = xunonz[colj+1]; 
  for(m=mstart; m<mstop; ++m) {
    if(rowu[m-1] == (rowi+1)) {
      unonz[m-1] += d;
      break;
    }
  }
}


template<class Scalar>
void
GenSpoolesSolver<Scalar>::add(GenFullM<Scalar> &kel, int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue;
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        // if(rowu[m-1] > unconstrNum[dofs[j]]+1)
        //   fprintf(stderr, "Bigger: %d %d\n", rowu[m-1]-1, unconstrNum[dofs[j]]);
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}


template<class Scalar>
void
GenSpoolesSolver<Scalar>::add(GenAssembledFullM<Scalar> &kel, int *dofs)
{
  // this function is used to assemble Kcc and requires dofs to be in constrained numbering
  int i, j, m, mstart, mstop, ri, rj;
  int kndof = kel.dim();                    // Dimension of element stiff.
  for(i = 0; i < kndof; ++i) {              // Loop over rows.
    if((ri = dofs[i]) == -1) continue;      // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {            // Loop over columns.
      if((rj = dofs[j]) == -1) continue;    // Skip constrained dofs
      if(rj < ri) continue;
      mstart = xunonz[rj];
      mstop  = xunonz[rj+1];
      for(m=mstart; m<mstop; ++m) {
        if(rowu[m-1] == (ri + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::add(GenFullM<Scalar> &knd, int fi, int fj)
{
  // PJSA 1-24-07 this needs to work for a rectangular matrix also
  int i, j, m, mstart, mstop, ri, rj;
  //int kndof = knd.dim();                           // Dimension of element stiff.
  for(i = 0; i < knd.numRow(); ++i ) {             // Loop over rows.
    if((ri =unconstrNum[fi+i]) == -1) continue;    // Skip constrained dofs
    for(j = 0; j < knd.numCol(); ++j) {            // Loop over columns.
      if((rj = unconstrNum[fj+j]) == -1) continue; // Skip constrained dofs
      if(rj < ri) continue;
      mstart = xunonz[rj];
      mstop  = xunonz[rj+1];
      for(m=mstart; m<mstop; ++m) {
        if(rowu[m-1] == (ri + 1) ) {
          unonz[m-1] += knd[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  // PJSA: 8-26-04 fixed this function, wasn't adding mass properly
  if(dof < 0) return;
  int cdof;
  if(unconstrNum) cdof = unconstrNum[dof]; // PJSA: dof is now in unconstrained numbering
  else cdof = dof;
  if(cdof < 0) return;

  int mstart = xunonz[cdof];
  int mstop  = xunonz[cdof+1];
  for(int m=mstart; m<mstop; ++m) {
    // if(rowu[m-1] > cdof+1)
    //   fprintf(stderr, "Bigger: %d %d\n", rowu[m-1]-1, cdof);
    if(rowu[m-1] == (cdof + 1)) {
      unonz[m-1] += dmass;
      break;
    }
  }

/*  OLD CODE  (this isn't correct)
  if(dof < 0) return;
  if(unconstrNum[dof] < 0) return;  // Skip constrained dofs
  int m = xunonz[unconstrNum[dof]];
  if(rowu[m-1] == (unconstrNum[dof] + 1)) {
    unonz[m-1] += dmass;
  }
*/
}

//PJSA warning THREADS defined causes bug for nonlinear
#define THREADS

template<class Scalar>
Scalar
GenSpoolesSolver<Scalar>::diag(int dof) const
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
/* PJSA this is done in symmetricScaling ... this function should be consistent with diag - 2
      if(unonz[m] == 0.0) {
        return (1.0);
      } else
*/
        return unonz[m];
    }
  }
  throw "GenSpoolesSolver<Scalar>::diag - 1 - this should never be reached";
}

template<class Scalar>
Scalar &
GenSpoolesSolver<Scalar>::diag(int dof)
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
        return unonz[m];
    }
  }
  throw "GenSpoolesSolver<Scalar>::diag - 2 - this should never be reached";
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::symmetricScaling()
{
  if(numUncon == 0) return;
  //HB: delete scale array if it already exists
  //    -> for NonLinSpooles, to avoid memory leak
  //    each time factor(...) is called (see reBuild(...))
  if(scale) { delete [] scale; scale = 0; }
  scale = new Scalar[numUncon];

  int i, j, m, mstart, mstop;
  for(i = 0; i < numUncon; ++i) {
    //scale[i] = Scalar(1.0) / ScalarTypes::sqrt(ScalarTypes::norm(diag(i)));
    Scalar d = ScalarTypes::norm(diag(i));
    scale[i] = (d != 0.0) ? Scalar(1.0)/ScalarTypes::sqrt(d) : Scalar(1.0); // PJSA 1-9-08 correction to allow for zero diagonals
  }
  for(j = 0; j < numUncon; ++j) {

    mstart = xunonz[j]-1;
    mstop  = xunonz[j+1]-1;
    for(m=mstart; m<mstop; ++m) {
      i = rowu[m] - 1;
      unonz[m] *= scale[i]*scale[j];
    }
  }
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::applyScaling(Scalar *vector)
{
  if(isScaled) {
    int i;
    for(i=0; i<numUncon; ++i)
      vector[i] *= scale[i];
  }
}


/****************************************************************************/
/****************************************************************************/
// Factor functions

template<class Scalar>
void
GenSpoolesSolver<Scalar>::factor()
{
  if(numUncon==0) return;
  //print();
  if(isScaled) symmetricScaling();
#ifdef THREADS
    numThreads = threadManager->numThr();
    allFactor(true);
#else
  numThreads = 1;
  allFactor(false); 
#endif
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::parallelFactor()
{
  if(numUncon==0) return;
  //print();
  if(isScaled) symmetricScaling();
  allFactor(true);
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::allFactor(bool fctIsParal)
{
#ifdef USE_SPOOLES
  int i, j;
  //double tt0 = getTime();
  // Scalar trace = 0;
  for(i = 0; i < neq; ++i)
  for(j = xunonz[i]; j <  xunonz[i+1]; ++j) {
    if(i >= rowu[j-1]-1) {
      InpMtx_inputRealOrComplexEntry(inpMtx, rowu[j-1]-1, i, unonz[j-1]);
      // trace += ( (rowu[j-1]%11)-5.5)* unonz[j-1];
    } else
      fprintf(stderr, "Weird\n");
  }

  //tt0 = getTime();
  InpMtx_changeStorageMode(inpMtx, INPMTX_BY_VECTORS);
  // OLD: InpMtx_changeStorageMode(inpMtx, INPMTX_BY_CHEVRONS);

/*
   -------------------------------------------------
   STEP 2 : find a low-fill ordering
   (1) create the Graph object
   (2) order the graph using multiple minimum degree
   -------------------------------------------------
*/

  graph = Graph_new();
  IVL *adjIVL = InpMtx_fullAdjacency(inpMtx);
  int nedges  = IVL_tsize(adjIVL);

  Graph_init2(graph, 0, neq, 0, nedges, neq, nedges, adjIVL, NULL, NULL);

// spooles_tau - upper bound on the magnitude of the largest element in L or U
// spooles_seed - random number seed used for ordering
// spooles_msglvl - message output level
// spooles_maxdomainsize - maximum subgraph size used by Spooles orderings. used to control the incomplete nested dissection 
//     process. Any subgraph whose weight is less than maxdomainsize is not split further
// spooles_maxzeros - maximum number of zeros allowed in a supernode/front
// spooles_maxsize - maximum number of internal columns in supernode/front

  int seed = domain->solInfo().spooles_seed;  // default is 532196
  int maxdomainsize = int(neq/domain->solInfo().spooles_maxdomainsize+0.5);  // default is neq/24
  int maxsize = domain->solInfo().spooles_maxsize;  // default is 64
  int maxzeros = int(neq*domain->solInfo().spooles_maxzeros+0.5);  // default is 0.04*neq
  //cerr << "neq = " << neq << ", maxdomainsize = " << maxdomainsize << ", maxsize = " << maxsize << ", maxzeros = " << maxzeros << ", seed = " << seed << endl;
  //cerr << "spooles_pivot = " << domain->solInfo().pivot << ", spooles_scale = " << domain->solInfo().spooles_scale << ", spooles_renum = " << domain->solInfo().spooles_renum << endl;

  //tt0=getTime();

  switch(domain->solInfo().spooles_renum) { // PJSA 10-11-2007
    default:
    case 0: // best of nested dissection and multisection ordering 
      frontETree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros, maxsize, seed, msglvl, msgfile); break;
    case 1: // multiple minimum degree
      frontETree = orderViaMMD(graph, seed, msglvl, msgfile); break;
    case 2: // multisection
      frontETree = orderViaMS(graph, maxdomainsize, seed, msglvl, msgfile); break;
    case 3: // nested dissection
      frontETree = orderViaND(graph, maxdomainsize, seed, msglvl, msgfile); break;
  }
  //frontETree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros, maxsize, seed, msglvl, msgfile); // best of nested dissection and multisection ordering

/*
   -----------------------------------------------------
   STEP 3: get the permutation, permute the matrix and
   -----------------------------------------------------
*/
  oldToNewIV = ETree_oldToNewVtxPerm(frontETree);
  int *oldToNew = IV_entries(oldToNewIV);
  newToOldIV = ETree_newToOldVtxPerm(frontETree);
  //int *newToOld = IV_entries(newToOldIV);

  ETree_permuteVertices(frontETree, oldToNewIV);
  InpMtx_permute(inpMtx, oldToNew, oldToNew);
  InpMtx_mapToUpperTriangle(inpMtx);
  InpMtx_changeCoordType(inpMtx, INPMTX_BY_CHEVRONS);
  InpMtx_changeStorageMode(inpMtx, INPMTX_BY_VECTORS);
  //tt0=getTime();

  symbfacIVL = SymbFac_initFromInpMtx(frontETree, inpMtx); // PJSA
  //tt0=getTime();

/*
   ------------------------------------------
   STEP 4: initialize the front matrix object
   ------------------------------------------
*/
  frontMtx   = FrontMtx_new();
  mtxManager = SubMtxManager_new();

#ifdef THREADS
  if(fctIsParal) {
    SubMtxManager_init(mtxManager, LOCK_IN_PROCESS, 0);
    FrontMtx_init(frontMtx, frontETree, symbfacIVL, SpoolesType<Scalar>::type, SPOOLES_SYMMETRIC,
                  FRONTMTX_DENSE_FRONTS, pivotingflag, LOCK_IN_PROCESS, 0, NULL,
                  mtxManager, msglvl, msgfile);
  } else
#endif
  {
    SubMtxManager_init(mtxManager, NO_LOCK, 0);
    FrontMtx_init(frontMtx, frontETree, symbfacIVL, SpoolesType<Scalar>::type, SPOOLES_SYMMETRIC,
                  FRONTMTX_DENSE_FRONTS, pivotingflag, NO_LOCK, 0, NULL,
                  mtxManager, msglvl, msgfile);
  }
/*
   -----------------------------------------
   STEP 5: compute the numeric factorization
   -----------------------------------------
*/
  ChvManager *chvmanager = ChvManager_new();
#ifdef THREADS
  if(fctIsParal)
    ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1);
  else
#endif
    ChvManager_init(chvmanager, NO_LOCK, 1);
  DVfill(22, cpus, 0.0); // PJSA changed 1st argument from 10 to 20
  IVfill(7, stats, 0);

  double tau = domain->solInfo().spooles_tau; // default is 100.0
  int error = 0;
  //double t0 = getTime();
  Chv *rootchv = NULL;
#ifdef THREADS
  if(fctIsParal) {
    int nfront;
    numThreads = threadManager->numThr();
    if(numThreads > (nfront=FrontMtx_nfront(frontMtx)))
      numThreads = nfront;
    cumopsDV = DV_new();
    DV_init(cumopsDV, numThreads, NULL);
    ownersIV = ETree_ddMap(frontETree, SpoolesType<Scalar>::type, SPOOLES_SYMMETRIC, cumopsDV,1./(2.*numThreads));
    rootchv = FrontMtx_MT_factorInpMtx(frontMtx, inpMtx, tau, 0.0, chvmanager,
                                           ownersIV, 0, &error, cpus, stats, msglvl, msgfile);
  } else
#endif
  rootchv = FrontMtx_factorInpMtx(frontMtx, inpMtx, tau, 0.0, chvmanager,
                                         &error, cpus, stats, msglvl, msgfile);
  if(rootchv != NULL) cerr << " ... WARNING: Matrix factored by spooles found to be singular ... \n";
  // note: spooles is only for full rank matrices. For symmetric positive semi definite matrices try sparse or skyline, and for general symmetric indefinite try mumps pivot

  ChvManager_free(chvmanager);

/*
   --------------------------------------
   STEP 6: post-process the factorization
   --------------------------------------
*/
  FrontMtx_postProcess(frontMtx, msglvl, msgfile) ;
  //int nf = FrontMtx_nfront(frontMtx);
  //SubMtx *sm = FrontMtx_diagMtx(frontMtx, nf-1);
  //FrontMtx_writeToFile(frontMtx,"frontMtx1");

#if DEBUG_SPOOLES >= 1
  cerr << " ... Spooles stats: \n"
       << "     # of pivots = " << stats[0] << endl
       << "     # of pivot tests = " << stats[1] << endl
       << "     # of delayed rows and cols = " << stats[2] << endl
       << "     # of entries in D = " << stats[3] << endl
       << "     # of entries in L = " << stats[4] << endl
       << "     # of entries in U = " << stats[5] << endl;
  cerr << " ... Spools timings: \n"
       << "     time to construct graph = " << cpus[0] << endl
       << "     time to compress graph = " << cpus[1] << endl
       << "     time to order graph = " << cpus[2] << endl
       << "     time for symbolic factorization = " << cpus[3] << endl
       << "     total setup time = " << cpus[4] << endl
       << "     time to setup the factorization = " << cpus[5] << endl
       << "     time to permute matrix = " << cpus[6] << endl
       << "     time to initialize front matrix = " << cpus[7] << endl
       << "     time to factor matrix = " << cpus[8] << endl
       << "     time to post-process matrix = " << cpus[9] << endl
       << "     total factor time = " << cpus[10] << endl;
#endif
#endif
}



/***********************************************************************************/
/***********************************************************************************/
// Solve functions

template<class Scalar>
void
GenSpoolesSolver<Scalar>::solve(Scalar *rhs)
{
  if(numUncon==0) return;
  Scalar *solution = (Scalar *) dbg_alloca(sizeof(Scalar)*numUncon);
  for(int i=0; i < numUncon; i++) solution[i] = 0.0;  // PJSA 8-25-04 debugging test
  solve(rhs, solution);
  for(int i=0; i < numUncon; i++)
    rhs[i] = solution[i];
}

/* PJSA: just for debugging
template<class Scalar>
double getMaxRealSq(Scalar *vec, int len) {
  double ret = 0.0;
  for(int i=0; i<len; ++i) {
    double tmp = ScalarTypes::Real(vec[i])*ScalarTypes::Real(vec[i]);
    if(tmp > ret) ret = tmp;
  }
  return ret;
}
*/

template<class Scalar>
void
GenSpoolesSolver<Scalar>::solve(Scalar *_rhs, Scalar *solution)
{
  this->solveTime = -getTime();
#ifdef USE_SPOOLES
  //double tt0 = getTime();
  Scalar *rhs = (Scalar *) dbg_alloca(sizeof(Scalar)*numUncon);
  for(int i=0; i < numUncon; i++)
    rhs[i] = _rhs[i];

  applyScaling(rhs);

  if(mtxB == 0) {
    mtxB = DenseMtx_new() ;
    DenseMtx_init(mtxB, SpoolesType<Scalar>::type, 0, 0, neq, 1, 1, neq);
  }
  DenseMtx_zero(mtxB) ;

  int i;
  for(i=0; i<neq; ++i)
    DenseMtx_setRealOrComplexEntry(mtxB, i, 0, rhs[i]);

/*
   ---------------------------------------------------------
   STEP 8: permute the right hand side into the new ordering
   ---------------------------------------------------------
*/
  DenseMtx_permuteRows(mtxB, oldToNewIV) ;

  if(mtxX == 0) {
    mtxX = DenseMtx_new();
    DenseMtx_init(mtxX, SpoolesType<Scalar>::type, 0, 0, neq, 1, 1, neq);
  }
  DenseMtx_zero(mtxX);

  //double tt1 = getTime();
#ifdef THREADS_SOLVE
  SolveMap *solvemap = SolveMap_new();
  SolveMap_ddMap(solvemap,SPOOLES_SYMMETRIC,FrontMtx_upperBlockIVL(frontMtx),
                 FrontMtx_lowerBlockIVL(frontMtx), numThreads, ownersIV, 
                 FrontMtx_frontTree(frontMtx), numThreads/2, 0, 0);
  FrontMtx_MT_solve(frontMtx, mtxX, mtxB, mtxManager, solvemap, cpus, msglvl, msgfile);
#else
  // FrontMtx_writeToFile(frontMtx,"frontMtx2");
  FrontMtx_solve(frontMtx, mtxX, mtxB, mtxManager, cpus, msglvl, msgfile); 
#endif

/*
   --------------------------------------------------------
   STEP 10: permute the solution into the original ordering
   --------------------------------------------------------
*/
  DenseMtx_permuteRows(mtxX, newToOldIV);

for(i=0;i<neq;i++) solution[i] = 0;

  // Retrieve the solution
  DenseMtx_realOrComplexEntries(mtxX, solution, neq);

  // scaling
  applyScaling(solution);

  // End Solution Phase

#if DEBUG_SPOOLES == 2
  cerr << " ... Spools timings: \n"
       << "     time to setup the parallel solve = " << cpus[0] << endl
       << "     time to permute rhs = " << cpus[1] << endl
       << "     time to solve = " << cpus[2] << endl
       << "     time to permute solution = " << cpus[3] << endl
       << "     total solve time = " << cpus[4] << endl;
#endif
#endif
  this->solveTime = +getTime();
}

template<class Scalar>
double
GenSpoolesSolver<Scalar>::getMemoryUsed(void)
{ 
  return 0;
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::print()
{
 cout << endl;
 cerr.setf(ios_base::scientific,ios_base::floatfield);
 cerr.precision(16);
 int i, mstart, mstop, m;
 for(i=0; i<numUncon; ++i) {
   mstart = xunonz[i];
   mstop  = xunonz[i+1];
   for(m=mstart; m<mstop; ++m)
     cerr << "K(" << i+1 << "," << rowu[m-1] << ") = " << unonz[m-1] << ",\n";
 }
}


template<class Scalar>
long
GenSpoolesSolver<Scalar>::size() 
{ 
  return neq;
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
  communicator->globalSum(xunonz[numUncon]-1, unonz);
#endif
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::zeroAll() 
{
#ifdef USE_SPOOLES
  cleanUp(); // PJSA 6-15-06
  // InpMtx_clearData(inpMtx);
  inpMtx = InpMtx_new();
  InpMtx_init(inpMtx, INPMTX_BY_COLUMNS, SpoolesType<Scalar>::type, nNonZero, neq);  // PJSA

  for(int i = 0; i < nNonZero; ++i)
    unonz[i] = 0.0;
#endif
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::cleanUp()
{
#ifdef USE_SPOOLES
  if(inpMtx) { InpMtx_free(inpMtx); inpMtx = 0; }
  if(frontMtx) { FrontMtx_free(frontMtx); frontMtx = 0; }
  if(newToOldIV) { IV_free(newToOldIV); newToOldIV = 0; }
  if(oldToNewIV) { IV_free(oldToNewIV); oldToNewIV = 0; }
  if(mtxManager) { SubMtxManager_free(mtxManager); mtxManager = 0; }
  if(ownersIV) { IV_free(ownersIV); ownersIV = 0; }
  if(mtxB) { DenseMtx_free(mtxB); mtxB = 0; }
  if(mtxX) { DenseMtx_free(mtxX); mtxX = 0; }
  if(symbfacIVL) { IVL_free(symbfacIVL); symbfacIVL = 0; } // PJSA
  if(frontETree) { ETree_free(frontETree); frontETree = 0; } // PJSA
  if(cumopsDV) { DV_free(cumopsDV); cumopsDV = 0; } // PJSA
  if(graph) { Graph_free(graph); graph = 0; } // PJSA
#endif
}

template<class Scalar>
GenSpoolesSolver<Scalar>::~GenSpoolesSolver()
{
  if(unonz) { delete [] unonz; unonz = 0; }
  if(scale) { delete [] scale; scale = 0; }
  if(msgfile) fclose(msgfile);
  cleanUp();
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::init()
{
  unonz = 0;
#ifdef USE_SPOOLES
  inpMtx = 0;
  frontMtx = 0;
  newToOldIV = 0;
  oldToNewIV = 0;
  mtxManager = 0;
  ownersIV = 0;
  mtxB = 0;
  mtxX = 0;
  symbfacIVL = 0;
  frontETree = 0;
  cumopsDV = 0; 
  graph = 0;
  pivotingflag = (domain->solInfo().pivot) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
#endif
  scale = 0;
  isScaled = domain->solInfo().spooles_scale;
  msglvl = domain->solInfo().spooles_msglvl;
  msgfile = (msglvl > 0) ? fopen("spooles_msgfile","w") : NULL;
  nrbm = 0;
  rbm = 0;
}


template<class Scalar>
void
GenSpoolesSolver<Scalar>::getRBMs(double *rigidBodyModes)
{
  if(rbm) rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::getRBMs(Vector *rigidBodyModes)
{
  if(rbm) rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSpoolesSolver<Scalar>::getRBMs(VectorSet& rigidBodyModes)
{
  if(rbm) rbm->getRBMs(rigidBodyModes);
}

template<>
void
GenSpoolesSolver<complex<double> >
   ::addImaginary(FullSquareMatrix &kel, int *dofs);

template<>
void
GenSpoolesSolver<double>
   ::addImaginary(FullSquareMatrix &kel, int *dofs);
