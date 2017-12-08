// ****************************************************************************************************
#include "SubDomain.h"
#include <Math.d/SparseSet.h>

template<class Scalar>
void
GenSubDomain<Scalar>::makeZstarAndR(double *centroid)
{
  rigidBodyModesG = new Rbm(dsa, c_dsa, nodes, sinfo.tolsvd,
                            centroid+3*group, cornerNodes, numCRN, numCRNdof, cornerDofs,
                            numMPC_primal, (SubLMPCons<Scalar> **) mpc_primal);
}

template<> 
void
GenSubDomain<double>::makeLocalRstar(FullM **Qtranspose) 
{
  FullM &R = rigidBodyModesG->R;
  FullM *Rc = rigidBodyModesG->Rc;
  if(numMPC_primal > 0) {
    FullM Qbody(Qtranspose[group]->transpose(), R.numCol(), bodyRBMoffset, Qtranspose[group]->numRow(), 0);
    Rstar = R * Qbody;
  }
  else {
    Rstar = R % *(Qtranspose[group]);
  }
}

template<>
void
GenSubDomain<DComplex>::makeLocalRstar(FullM **Qtranspose)
{
  std::cerr << "GenSubDomain<DComplex>::makeLocalRstar(FullM **Qtranspose) is not implemented\n";
}

template<class Scalar>
void
GenSubDomain<Scalar>::useKrrNullspace()
{
  // EXPERMENTAL... use alaebraic null space of Krr with no corners
  int neq = this->Krr->neqs();
  int nzem = this->Krr->numRBM();
  Rstar.setNewSize(neq, nzem);
  if(nzem > 0) {
    std::vector<Scalar> rbmv(neq*nzem);
    this->Krr->getNullSpace(rbmv.data());

    // Copy rigid body modes (rbm) into Rstar
    for(int m=0; m<nzem; ++m) {
      for(int i=0; i<neq; ++i)
        Rstar[i][m] = rbmv[i+m*neq];
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::addRalpha(Scalar *u, GenVector<Scalar> &alpha) const
{
  int i, j;
  for(i=0; i<Rstar.numRow(); ++i)
    for(j=0; j<Rstar.numCol(); ++j)
      u[i] += Rstar[i][j] * alpha[groupRBMoffset + j];
}

template<class Scalar>
void
GenSubDomain<Scalar>::addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur) const
{
  int numCDofs = (Src) ? Src->numCol() : 0;
  Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
  Scalar *localr = new Scalar[localLen()];

  for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
    Scalar *rbm = rbms + glNumCDofs*iRbm;
    for(int i=0; i<numCDofs; ++i) {
      if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]*alpha[iRbm]; }
      else localc[i] = 0.0;
    }

    for(int i = 0; i < localLen(); ++i) localr[i] = 0.0;
    if(Src) Src->transposeMultAdd(localc, localr);
    if(this->Krr) this->Krr->reSolve(localr);

    for(int i = 0; i < localLen(); ++i) ur[i] -= localr[i];
  }

  delete [] localr;
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleE(GenVector<Scalar> &e, Scalar *f) const
{
  if(numGroupRBM > 0) {
    GenVector<Scalar> local_e(numGroupRBM, 0.0);
    GenVector<Scalar> fvec(f, Rstar.numRow());
    local_e = Rstar ^ fvec; // = Rtranspose * fvec
    e.add(local_e, groupRBMoffset);  
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleTrbmE(Scalar *rbms, int nrbms, int size, Scalar *e, Scalar *fr) const
{
  int numCDofs = (Src) ? Src->numCol() : 0;
  Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
  for(int i=0; i<numCDofs; ++i) localc[i] = 0.0;

  Scalar *localr = new Scalar[localLen()];
  for(int i = 0; i < localLen(); ++i) localr[i] = -fr[i];
  if(this->Krr) this->Krr->reSolve(localr);
  if(Src) Src->multAdd(localr, localc); // localc = - (Krr^-1 Krc)^T fr
  delete [] localr;

  for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
    Scalar *rbm = rbms + size*iRbm;
    for(int i=0; i<numCDofs; ++i) {
      if(cornerEqNums[i] > -1) { e[iRbm] += rbm[cornerEqNums[i]]*localc[i]; } // e += -N^T (Krr^-1 Krc)^T fr
    }
  }
}

// ****************************************************************************************************

template<class Scalar>
void
GenSubDomain<Scalar>::makeG()
{
  // make G for each potential contact/mpc neighbour
  G = new GenFullM<Scalar> * [scomm->numT(SComm::mpc)];
  neighbG = new GenFullM<Scalar> * [scomm->numT(SComm::mpc)];
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {  // loop through all potential contact/mpc neighbours
    neighbG[i] = 0; 
    G[i] = new GenFullM<Scalar>(scomm->lenT(SComm::mpc,i), numGroupRBM);
    if(numGroupRBM > 0) {
      G[i]->zero();
      for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {  // loop through potential contact/mpc nodes
        int locMpcNb = scomm->mpcNb(i,j);
        SubLMPCons<Scalar> *m = mpc[locMpcNb];
        for(int k = 0; k < m->nterms; k++) {
          int cDof = (m->terms)[k].cdof;
          if(cDof > -1) {
            for(int iRbm = 0; iRbm < numGroupRBM; ++iRbm)
              (*G[i])[j][iRbm] += Rstar[cDof][iRbm]*(m->terms)[k].coef;
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeTrbmG(Scalar *rbms, int nrbm, int size)
{
  // rbms is the null space of the global Kcc^* matrix
  // nrbm is the nullity of the global Kcc^* matrix
  // size is the number of rows and columns of the global Kcc^* matrix
  // TODO what about augmentation
  int numCDofs = (Src) ? Src->numCol() : 0;
  Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);

  int nrbms_local = 0; int first = 0;
  std::map<int, int> localToGlobalRBM;
  std::map<int, int> globalToLocalRBM;
  for(int iRbm = 0; iRbm < nrbm; ++iRbm) {
    Scalar *rbm = rbms + size*iRbm;
    Scalar dot = 0.0;
    for(int i=0; i<numCDofs; ++i) {
      if(cornerEqNums[i] > -1) { dot += rbm[cornerEqNums[i]]*rbm[cornerEqNums[i]]; }
    }
    if(dot != 0.0) {
      localToGlobalRBM[nrbms_local] = iRbm;
      globalToLocalRBM[iRbm] = nrbms_local;
      if(nrbms_local == 0) first = iRbm;
      nrbms_local++;
    }
  }

  numGroupRBM = nrbms_local; // TODO this isn't general since global trbms may not be grouped like grbms
  groupRBMoffset = first;    // TODO this isn't general

  G = new GenFullM<Scalar>*[scomm->numT(SComm::mpc)];
  neighbG = new GenFullM<Scalar>*[scomm->numT(SComm::mpc)];
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    neighbG[i] = 0;
    G[i] = new GenFullM<Scalar>(scomm->lenT(SComm::mpc,i), nrbms_local);
    G[i]->zero();
  }

  if(nrbms_local == 0 || numMPC == 0) return;

  Scalar *localr = new Scalar[localLen()];

  for(int iRbm = 0; iRbm < nrbms_local; ++iRbm) {
    int glRbmId = localToGlobalRBM[iRbm];
    Scalar *rbm = rbms + size*glRbmId;
    for(int i=0; i<numCDofs; ++i) {
      if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]; }
      else localc[i] = 0.0;
    }

    // G = (-Br^(s)Krr^{-1}Krc + Bc)Lcc Nc
    for(int i=0; i<localLen(); ++i) localr[i] = 0.0;
    if(Src) Src->transposeMultSubtract(localc, localr);
    if(this->Krr) this->Krr->reSolve(localr);
    for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
      for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
        int locMpcNb = scomm->mpcNb(i,j);
        SubLMPCons<Scalar> *m = mpc[locMpcNb];
        for(int k = 0; k < m->nterms; k++) {
          int cc_dof = (m->terms)[k].ccdof;
          if(cc_dof >= 0) (*G[i])[j][iRbm] += localr[cc_dof]*(m->terms)[k].coef;
          else {
            int dof = (m->terms)[k].dof;
            if((dof >= 0) && (cornerMap[dof] >= 0))
              (*G[i])[j][iRbm] += localc[cornerMap[dof]] * (m->terms)[k].coef;
          }
        }
      }
    }
  }
  delete [] localr;
}

template<class Scalar>
void
GenSubDomain<Scalar>::setGCommSize(FSCommStructure *pat) const
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int nRow = G[i]->numRow();
    int nCol = numGroupRBM;
    pat->setLen(subNumber, scomm->neighbT(SComm::mpc, i), nRow*nCol);
  }   
}
 
template<class Scalar>
void
GenSubDomain<Scalar>::sendG(FSCommPattern<Scalar> *rbmPat)
{
  if(numGroupRBM == 0) return;
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) 
    rbmPat->sendData(subNumber, scomm->neighbT(SComm::mpc, i),  G[i]->data());
}

template<class Scalar>
void
GenSubDomain<Scalar>::receiveG(FSCommPattern<Scalar> *rbmPat)
{
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    FSSubRecInfo<Scalar> rInfo = rbmPat->recData(scomm->neighbT(SComm::mpc, i), subNumber);
    int nRow = G[i]->numRow();  // number of potential contact dofs on interface with neighb
    int nCol = neighbNumGroupGrbm[i];  //number of rbms for neighb's group
    if(neighbG[i]) delete neighbG[i];
    neighbG[i] = new GenFullM<Scalar>(rInfo.data, nRow, nCol);  
  }
}

template<class Scalar>
void GenSubDomain<Scalar>::zeroG()
{
  if(G) {
    for(int i=0; i<scomm->numT(SComm::mpc); ++i)
      if(G[i]) G[i]->zero();
  }
  if(neighbG) {
    for(int i=0; i<scomm->numT(SComm::mpc); ++i)
      if(neighbG[i]) neighbG[i]->zero();
  }
}

template<class Scalar>
void GenSubDomain<Scalar>::deleteG()
{
  if(G) {
    for(int i=0; i<scomm->numT(SComm::mpc); ++i)
      if(G[i]) { delete G[i]; G[i] = 0; }
     delete [] G; G = 0;
  }
  if(neighbG) {
    for(int i=0; i<scomm->numT(SComm::mpc); ++i)
      if(neighbG[i]) { delete neighbG[i]; neighbG[i] = 0; }
     delete [] neighbG; neighbG = 0;
  }
  if(neighbNumGroupGrbm) { delete [] neighbNumGroupGrbm; neighbNumGroupGrbm = 0; }
  if(neighbGroupGrbmOffset) { delete [] neighbGroupGrbmOffset; neighbGroupGrbmOffset = 0; }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multG(const GenVector<Scalar> &x, Scalar *y, Scalar alpha) const
{
  // y += alpha * G * x
  Scalar *mpcvec = new Scalar[numMPC];
  for(int i = 0; i < numMPC; ++i) mpcvec[i] = 0.0;
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int neighb = scomm->neighbT(SComm::mpc, i);
    for(int j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
      int locMpcNb = scomm->mpcNb(i,j);
      if(mpcvec[locMpcNb] == 0.0)
        for(int k = 0; k < numGroupRBM; ++k)
          mpcvec[locMpcNb] += (*G[i])[j][k] * x[k + groupRBMoffset];
      if(subNumber != neighb)
        for(int k = 0; k < neighbNumGroupGrbm[i]; ++k)
          mpcvec[locMpcNb] += (*neighbG[i])[j][k] * x[k + neighbGroupGrbmOffset[i]];
    }
  }
  for(int i = 0; i < scomm->lenT(SComm::mpc); ++i)
    y[scomm->mapT(SComm::mpc,i)] += alpha*mpcvec[scomm->boundDofT(SComm::mpc,i)];
  delete [] mpcvec;
}

template<class Scalar>
void 
GenSubDomain<Scalar>::trMultG(const Scalar *x, GenVector<Scalar> &y, Scalar alpha) const
{
  // compute y += alpha * G^t * x
  bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
      int locMpcNb = scomm->mpcNb(i,j);
      int iDof = scomm->mapT(SComm::mpc,i,j);
      if(mpcFlag[locMpcNb]) {
        for(int k = 0; k < numGroupRBM; ++k)
          y[k + groupRBMoffset] += alpha * (*G[i])[j][k] * x[iDof];
        mpcFlag[locMpcNb] = false;
      }
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver)
{
  if(numGroupRBM == 0) return;
  bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
  for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
    int numGroupRBM2 = neighbNumGroupGrbm[i];
    int groupRBMoffset2 = neighbGroupGrbmOffset[i];
    GenVector<Scalar> d(scomm->lenT(SComm::mpc,i));
    for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
      int locMpcNb = scomm->mpcNb(i,j);
      d[j] = (mpc[locMpcNb]->active) ? 0.0 : 1.0;
    }
    if((numGroupRBM2 > 0) && (subNumber != scomm->neighbT(SComm::mpc,i))) {
      GenFullM<Scalar> tmp2(numGroupRBM, numGroupRBM2);  // coupling term
      G[i]->transposeMultD(*neighbG[i], d, tmp2); // tmp2 = G^T * D * neighbG
      GtGsolver->add(tmp2, groupRBMoffset, groupRBMoffset2);
    }
    for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
      int locMpcNb = scomm->mpcNb(i,j);
      if(!mpcFlag[locMpcNb]) d[j] = 0.0; // prevents duplication for mpc shared between more than 2 subdomains
      else mpcFlag[locMpcNb] = false;
    }
    GenFullM<Scalar> tmp(numGroupRBM, numGroupRBM);
    G[i]->transposeMultD(*G[i], d, tmp); // tmp = G^T * D * G
    GtGsolver->add(tmp, groupRBMoffset, groupRBMoffset);
  }
}

// ************************************************************************************************

template<class Scalar>
void
GenSubDomain<Scalar>::buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, const Connectivity *cornerToSub)
{
  int i,j,k;
  if(numGroupRBM == 0) numGlobalRBMs = 0;
  else {
    numGlobalRBMs = Xmatrix.numCol();
    GenFullM<Scalar> groupX(Xmatrix, numGroupRBM, groupRBMoffset, numGlobalRBMs, 0); 
    Rstar_g = Rstar * groupX;

    if(!cornerToSub) return;
    double *sharedUse = new double[Rstar.numRow()];
    for(i=0; i<Rstar.numRow(); ++i) sharedUse[i] = 1.0;
    // if i is a shared dof set sharedUse[i] = 0.0
    // for all but one of the subdomains sharing it in this body

    for(i=0; i<scomm->numT(SComm::std); ++i) {  // check non-corner dofs
      if(subNumber > scomm->neighbT(SComm::std,i))
        for(j=0; j<scomm->lenT(SComm::std,i); ++j) 
          sharedUse[ccToC[scomm->boundDofT(SComm::std,i,j)]] = 0.0;
    }
    for(i=0; i<numCRN; ++i) { // check corner dofs
      if(subNumber != (*cornerToSub)[glCornerNodes[i]][0]) {
        int lDof[6];
        c_dsa->number(cornerNodes[i], cornerDofs[i].list(), lDof);
        for(j=0; j<cornerDofs[i].count(); ++j)
          if(lDof[j] >= 0) sharedUse[lDof[j]] = 0.0;
      }
    }

    // if i is a shared dof sharedRstar_g[i] is set to zero 
    // for all but one of the subdomains sharing it in this body
    // (used to prevent duplication in construction of RtR)
    if(sharedRstar_g) delete sharedRstar_g;
    sharedRstar_g = new GenFullM<Scalar>(Rstar_g.numRow(), Rstar_g.numCol());
    for(i=0; i<Rstar_g.numRow(); ++i) 
      for(j=0; j<Rstar_g.numCol(); ++j) (*sharedRstar_g)[i][j] = Rstar_g[i][j] * sharedUse[i];
    delete [] sharedUse;
  
    // if i is a shared dof tmpRstar_g[i] is set to Rstar_g[i]/n 
    // where n is the number of subdomains (is this body) sharing this dof
    // used to compute a distributed vector, since distvec[i]*n = actual value at dof i  
    if(tmpRstar_g) delete tmpRstar_g;
    tmpRstar_g = new GenFullM<Scalar>(Rstar_g);
    for(i=0; i<scomm->numT(SComm::all); ++i) {
      for(j=0; j<scomm->lenT(SComm::all,i); ++j) {
        int bdof = scomm->boundDofT(SComm::all,i,j);
        if(bdof >= 0)  // not a contact dof 
          for(k=0; k<numGlobalRBMs; ++k) 
            (*tmpRstar_g)[ccToC[bdof]][k] /= weight[bdof]; 
      }
    }
    for(i=0; i<numCRN; ++i) {
      int lDof[6];
      c_dsa->number(cornerNodes[i], cornerDofs[i].list(), lDof);
      for(j=0; j<cornerDofs[i].count(); ++j)
        if(lDof[j] >= 0)
          for(k=0; k<numGlobalRBMs; ++k) 
            (*tmpRstar_g)[lDof[j]][k] /= (double) (cornerToSub->num(glCornerNodes[i]));
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::getGlobalRBM(int iRBM, Scalar *Rvec) const
{
  if(numGlobalRBMs > 0) 
    for(int iRow=0; iRow<Rstar_g.numRow(); ++iRow)
	    Rvec[iRow] = Rstar_g[iRow][iRBM];
}

template<class Scalar>
void
GenSubDomain<Scalar>::subtractRstar_g(Scalar *u, GenVector<Scalar> &beta) const
{
	int i;
	if(numGlobalRBMs > 0) {
		// compute u = u - Rstar_g * beta  (second part of displacement projection)
		GenVector<Scalar> tmpu(Rstar_g.numRow());
		tmpu = Rstar_g * beta;
		for(i=0; i<Rstar_g.numRow(); ++i)
			u[i] = u[i] - tmpu[i];
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::addRstar_gT(Scalar *u, GenVector<Scalar> &beta) const
{
  if(numGlobalRBMs > 0) {
    // compute beta += Rstar_g^t * u  (first part of displacement projection)
    GenVector<Scalar> uvec(u, Rstar_g.numRow());
    beta += *tmpRstar_g ^ uvec;
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleRtR(GenFullM<Scalar> &RtR)
{
  // builds RtR for u projection
  if(numGlobalRBMs > 0) {
    GenFullM<Scalar> tmp(numGlobalRBMs, numGlobalRBMs);
    sharedRstar_g->transposeMult((*sharedRstar_g), tmp);
    RtR.add(tmp, 0, 0);
  }
}

// ****************************************************************************************************