#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Math.d/DiagMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Utils.d/Memory.h>
#include <Interface.d/MpcLocal.h>
#include <Element.d/Element.h>

template<class Scalar>
GenSandiaSubD<Scalar>::GenSandiaSubD(int _globalSubNumber, int _localSubNumber, 
                                     int _numnodes, const double *x, const double *y,
                                     const double *z, int *_glNums, unsigned short int *_dofMap,
                                     int _numele, int *eptr, int *econ, const int *_map, int *_subMap) 
 : GenSubDomain<Scalar>(_globalSubNumber, _localSubNumber)
{
  initialize();

  this->numnodes = _numnodes;
  this->numele = _numele;
  this->glNums = _glNums;
  dofMap = _dofMap;
  map = _map;
  subMap = _subMap;

  int i;
  for(i = 0; i < this->numnodes; ++i) {
    double xyz[3];
    xyz[0] = x[i];
    xyz[1] = y[i];
    xyz[2] = z[i];
    this->nodes.nodeadd(i, xyz);
  }
  this->makeGlobalToLocalNodeMap();
 
  for(i = 0; i < this->numele; ++i) {
    Element *ele= new SandiaElem(eptr[i+1]- eptr[i], econ+eptr[i], dofMap);
    this->packedEset.elemadd(i, ele);
  } 

  this->sinfo = domain->solInfo();
  this->sinfo.renum = 1; // Sloan renumbering
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::initialize()
{
  memoryK = 0; memoryPrec = 0; 
  this->numMPC   = 0; this->mpc = 0;
  rbm = 0;
  subSysMatrix = false;
  sandiaToLocalMap = 0;
  numContactNeighb = 0; contactNeighb = 0; neighbPtr = 0;
  mpcNeighb = 0;
}

template<class Scalar>
GenSandiaSubD<Scalar>::~GenSandiaSubD()
{
  if(sandiaToLocalMap) { delete [] sandiaToLocalMap; sandiaToLocalMap = 0; }
  if(subMap) { delete [] subMap; subMap = 0; }
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::initSysMatrix(int neq, GenSolver<Scalar> *&smat, GenSparseMatrix<Scalar> *&diagmat)
{
  if(subSysMatrix) deleteSysMatrix(smat, diagmat);

  if(this->numMPC) this->makeKbbMpc();
  else this->makeKbb(this->getCCDSA());

  // create the map from 'natural indexing' to our dsa numbering
  int *myMap = new int[neq];
  int mp = 0;
  int iNode, iDof;
  if(subMap) for(int i = 0; i < neq; ++i) myMap[i] = -1;
  for(iNode = 0; iNode < this->numnodes; ++iNode) {
    int dof = this->dsa->firstdof(iNode);
    int ndof = this->dsa->weight(iNode);
    for(iDof = 0; iDof < ndof; ++iDof) {
      if(this->c_dsa->getRCN(dof) >= 0) {
        if(subMap) myMap[subMap[mp++]] = dof;
        else myMap[mp++] = dof;
      }
      dof++;
    }
  }

  // Now create the map from Incoming indexing to dsa numbering
  sandiaToLocalMap = new int[neq];
  for(iDof = 0; iDof < neq; ++iDof) {
    sandiaToLocalMap[iDof] = myMap[map[iDof]];
  }
  delete [] myMap;

  diagmat = new GenDiagMatrix<Scalar>(this->cc_dsa);

  this->zeroEdgeDofSize(); // PJSA: for rebuilding Q
  this->makeQ();
  this->constructKcc();
  this->constructKrc();

  // ... Construct geometric rigid body modes if necessary ??

  dbg_alloca(0);
  switch(this->solInfo().getFetiInfo().solvertype) {
    default:
    case FetiInfo::sparse:
      smat = this->template constructBLKSparseMatrix<Scalar>(this->cc_dsa);
      smat->zeroAll();
      break;
    case FetiInfo::skyline:
      smat = this->template constructSkyMatrix<Scalar>(this->cc_dsa);
      break;
    case FetiInfo::spooles:
      smat = this->template constructSpooles<Scalar>();
      break;
    case FetiInfo::mumps:
      smat = this->template constructMumps<Scalar>();
      break;
  }

  if(this->c_dsa->size() > 0 && (this->c_dsa->size() - this->dsa->size()) != 0)
    this->Kuc = new GenCuCSparse<Scalar>(this->nodeToNode, this->dsa, this->c_dsa);
  else
    this->Kuc = 0;

  if(this->Krc) this->Krc->zeroAll();
  subSysMatrix = true;
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::deleteSysMatrix(GenSolver<Scalar> *smat, GenSparseMatrix<Scalar> *diagmat)
{
  // deallocate everything that is allocated in initSubSysMatrix(...)
  // not Krc !!!
  if(diagmat) { delete diagmat; diagmat = 0; }
  if(smat) { delete smat; smat = 0; this->KrrSparse = 0; this->Krr = 0; }
  if(this->KiiSolver) { delete this->KiiSolver; this->KiiSolver = 0; this->KiiSparse = 0; }
  if(this->Kib) { delete this->Kib; this->Kib = 0; }
  if(this->Kbb) { delete this->Kbb; this->Kbb = 0; }
  if(this->Kcc) { delete this->Kcc; this->Kcc = 0; }
  if(this->Kuc) { delete this->Kuc; this->Kuc = 0; }
  if(this->fcstar) { delete [] this->fcstar; this->fcstar = 0; }
  if(sandiaToLocalMap) { delete [] sandiaToLocalMap; sandiaToLocalMap = 0; }
  subSysMatrix = false;
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::zeroSysMatrix(GenSolver<Scalar> *smat, GenSparseMatrix<Scalar> *diagmat)
{
  if(diagmat) ((GenDiagMatrix<Scalar>*) diagmat)->zeroAll();
  if(this->KiiSparse) this->KiiSparse->zeroAll();
  if(this->Kib) this->Kib->zeroAll();
  if(this->Kbb) this->Kbb->zeroAll();
  if(this->Krc) this->Krc->zeroAll();
  if(this->Kcc) {
    this->Kcc->zero();
    for(int i=0; i<this->Kcc->dim(); ++i) this->fcstar[i] = 0.0;
  }
  if(smat) smat->zeroAll();
  if(this->Kuc) this->Kuc->zeroAll();
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::addSysMatrix(int neq, const double *Kaa, const int *Kaj,
                                    const int *Kai, Scalar multiplier, 
                                    GenSolver<Scalar> *&smat, GenSparseMatrix<Scalar> *&diagmat, bool isK)
{
  if(isK || (this->sinfo.getFetiInfo().prectype == FetiInfo::shifted)) { // assemble preconditioner
    if(this->KiiSolver) this->KiiSolver->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
    if(this->Kib) this->Kib->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
    if(this->Kbb) this->Kbb->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
  }
  //else cerr << "Preconditioner is not shifted\n";
  
  if(diagmat) ((GenDiagMatrix<Scalar>*) diagmat)->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
  if(smat) smat->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
  if(this->Krc) this->Krc->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
  if(this->Kcc) this->Kcc->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
  if(this->Kuc) this->Kuc->addBoeing(neq, Kai, Kaj, Kaa, sandiaToLocalMap, multiplier);
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::makeF(int neq, Scalar *f, Scalar *lf)
{
  // salinas to feti mapping
  for(int iDof = 0; iDof < neq; ++iDof) {
    if(sandiaToLocalMap[iDof] > -1)
      lf[this->c_dsa->getRCN(sandiaToLocalMap[iDof])] = f[iDof];
  }
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::getD(int neq, Scalar *dis, Scalar *ldis, bool assemble)
{
  // feti to salinas mapping
  if(assemble) {
    for(int iDof = 0; iDof < neq; ++iDof) {
      if(sandiaToLocalMap[iDof] > -1)
        dis[iDof] += ldis[this->c_dsa->getRCN(sandiaToLocalMap[iDof])];
    }
  }
  else {
    for(int iDof = 0; iDof < neq; ++iDof) {
      if(sandiaToLocalMap[iDof] > -1)
        dis[iDof] = ldis[this->c_dsa->getRCN(sandiaToLocalMap[iDof])];
    }
  }
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::remap(int neq, Scalar *salinasDis, Scalar *fetiDis)
{
  for(int iDof = 0; iDof < neq; ++iDof) {
    if(sandiaToLocalMap[iDof] > -1)
      salinasDis[iDof] = fetiDis[this->c_dsa->getRCN(sandiaToLocalMap[iDof])];
  }
}

// HB: from salinas dof ordering to Feti dof ordering (inverse of GenSandiaSubD<Scalar>::remap)
// same as GenSandiaSubD<Scalar>::makeF ... should make only one to be used for both ...
// Should replace double* by Scalar* in  near future ...
template<class Scalar>
double
GenSandiaSubD<Scalar>::salinasToFetiVector(double *salinasDisp, double *fetiDis)
{
  int neq = this->c_dsa->size();
  double norm = 0.0;
  for(int iDof = 0; iDof < neq; ++iDof) fetiDis[iDof] = 0.0; // maybe not worthy ...
  for(int iDof = 0; iDof < neq; ++iDof) {
    fetiDis[this->c_dsa->getRCN(sandiaToLocalMap[iDof])] = salinasDisp[iDof];
    norm += salinasDisp[iDof]*salinasDisp[iDof];
  }
  return(sqrt(norm));
}


template<class Scalar>
double
GenSandiaSubD<Scalar>::fetiToSalinasVector(double* fetiDis, double* salinasDis)
{
  int neq = this->c_dsa->size();
  for(int iDof = 0; iDof < neq; ++iDof)
    salinasDis[iDof] = fetiDis[this->c_dsa->getRCN(sandiaToLocalMap[iDof])];
  return(0.0);
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::getRBMs(int size, int nRBM, Scalar *salinasRBMs, Scalar *fetiRBMs)
{
  int neq = this->c_dsa->size();
  for(int iMode = 0; iMode < nRBM; ++iMode) {
    Scalar *salinasRBM = salinasRBMs + iMode*size;
    Scalar *fetiRBM = fetiRBMs + iMode*size;
    remap(neq, salinasRBM, fetiRBM);
  }
}

template<class Scalar>
void
GenSandiaSubD<Scalar>::setMpcGlobalTermIds(MpcLocal *mpclocal, int *global_node_nums, bool oneSubPerCPU)
{
  for(int iMPC=0; iMPC<this->numMPC; ++iMPC) {
    int glMpcNb = this->localToGlobalMPC[iMPC];
    this->mpc[iMPC]->gsize = mpclocal[glMpcNb].NumEntriesGlobal();
    for(int i=0; i<this->mpc[iMPC]->nterms; ++i) {
       if(oneSubPerCPU) {
         this->mpc[iMPC]->gi[i] = mpclocal[glMpcNb].LocalToGlobalEntry(i);
       }
       else { // mpclocal may be further split between multiple subdomains on this cpu
         for(int j=0; j<mpclocal[glMpcNb].NumEntries(); ++j) { 
           if(global_node_nums[mpclocal[glMpcNb].LocalId(j)] == this->glNums[this->mpc[iMPC]->terms[i].nnum] && 
              mpclocal[glMpcNb].NodalDof(j) == this->mpc[iMPC]->terms[i].dofnum) {
              this->mpc[iMPC]->gi[i] = mpclocal[glMpcNb].LocalToGlobalEntry(j);
              break;
           }
         }
       }
    }
  }
}
