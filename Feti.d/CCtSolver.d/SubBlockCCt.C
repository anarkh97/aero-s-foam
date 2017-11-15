template<class Scalar>
SubBlockCCtSolver<Scalar>::SubBlockCCtSolver(Connectivity *mpcToMpc, Connectivity *_mpcToSub, 
                                             int _numSubsWithMpcs, GenSubDomain<Scalar> **_subsWithMpcs,
                                             FSCommunicator *_fetiCom, Connectivity *_cpuToSub)
{
  mpcToSub = _mpcToSub;
  this->numSubsWithMpcs = _numSubsWithMpcs;
  this->subsWithMpcs = _subsWithMpcs;
  this->fetiCom = _fetiCom;
  cpuToSub = _cpuToSub;
  myCPU = _fetiCom->cpuNum();
  
  paralApply(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::makeLocalMpcToGlobalMpc, mpcToMpc);
  mpcvPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSubsWithMpcs; ++i) this->subsWithMpcs[i]->setMpcCommSize(mpcvPat);  // this is used for combineMpcInterfaceVec()
  mpcvPat->finalize();
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::constructLocalCCtsolver);
  cctPat = 0;
}

template<class Scalar>
SubBlockCCtSolver<Scalar>::~SubBlockCCtSolver()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::deleteLocalCCtsolver);
  delete mpcvPat;
  delete cctPat;
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::assemble()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::assembleLocalCCtsolver);
  if(!cctPat) {
    cctPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                       FSCommPattern<Scalar>::NonSym);
    for(int i=0; i<this->numSubsWithMpcs; ++i) this->subsWithMpcs[i]->setCCtCommSize(cctPat);
    cctPat->finalize();
  }
  paralApply(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::sendNeighbCCtsolver, cctPat, mpcToSub);
  cctPat->exchange();
  paralApply(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::recNeighbCCtsolver, cctPat, mpcToSub);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::factor()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::factorLocalCCtsolver);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::zeroAll()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &GenSubDomain<Scalar>::zeroLocalCCtsolver);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::reSolve(GenDistrVector<Scalar> &v)
{
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::solveLocalCCt, v);
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::sendMpcInterfaceVec, v);
  mpcvPat->exchange();
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::combineMpcInterfaceVec, v);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::solveLocalCCt(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->solveLocalCCt(subv);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::sendMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->sendMpcInterfaceVec(mpcvPat, subv);
}
template<class Scalar>
void
SubBlockCCtSolver<Scalar>::combineMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->combineMpcInterfaceVec(mpcvPat, subv);
}
