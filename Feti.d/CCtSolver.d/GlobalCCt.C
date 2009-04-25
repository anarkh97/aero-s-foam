
extern Domain * domain;

template<class Scalar>
GlobalCCtSolver<Scalar>::GlobalCCtSolver(Connectivity *mpcToMpc, Connectivity *_mpcToCpu, int _numSubsWithMpcs, 
                                         GenSubDomain<Scalar> **_subsWithMpcs, FetiInfo *finfo, FSCommunicator *_fetiCom)
{
  filePrint(stderr, " ... Build CCt solver               ...\n");
  this->mpcToCpu = _mpcToCpu;
  this->numSubsWithMpcs = _numSubsWithMpcs;
  this->subsWithMpcs = _subsWithMpcs;
  this->fetiCom = _fetiCom;
  this->glNumMpc = mpcToMpc->csize();
  switch(finfo->cctSolver) {
    default:
    case (FetiInfo::sparse) : 
      {
        mpcEqNums = new SimpleNumberer(this->glNumMpc);
        for(int i = 0; i < this->glNumMpc; ++i) mpcEqNums->setWeight(i, 1);
        mpcEqNums->makeOffset();
        CCtsolver = new GenBLKSparseMatrix<Scalar>(mpcToMpc, mpcEqNums, finfo->cct_tol, domain->solInfo().sparse_renum);
        CCtsolver->zeroAll();
      }
      break;
    case (FetiInfo::skyline) :
      {
        compStruct renumber = mpcToMpc->renumByComponent(1); // sloan renumbering
        //if(this->fetiCom->cpuNum() == 0) cerr << "CCt has " << renumber.numComp << " blocks\n";
        mpcEqNums = new SimpleNumberer(this->glNumMpc, renumber.renum, 1);
        for(int i = 0; i < this->glNumMpc; ++i) mpcEqNums->setWeight(i, 1);
        mpcEqNums->makeOffset();
        int scaledFlag = (finfo->cctScaled) ? 1 : 0;
        CCtsolver = new GenSkyMatrix<Scalar>(mpcToMpc, mpcEqNums, finfo->cct_tol, scaledFlag);
        delete [] renumber.xcomp;
      }
      break;
#ifdef USE_SPOOLES
      case FetiInfo::spooles: {
        mpcEqNums = new SimpleNumberer(this->glNumMpc);
        for(int i = 0; i < this->glNumMpc; ++i) mpcEqNums->setWeight(i, 1);
        mpcEqNums->makeOffset();
        CCtsolver = new GenSpoolesSolver<Scalar>(mpcToMpc, mpcEqNums);
      } break;
#endif

  }
  CCtsolver->setPrintNumTrbm(false);
}

template<class Scalar>
GlobalCCtSolver<Scalar>::~GlobalCCtSolver()
{
  delete mpcEqNums;
  delete CCtsolver;
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::assemble()
{
  bool new_global_cct = true;
  if(new_global_cct && (this->numSubsWithMpcs > 1)) { 
     // each sub (with LMPCs) compute its own contributions to the global CCt matrix
     execParal1R(this->numSubsWithMpcs, this, &GlobalCCtSolver<Scalar>::computeSubContributionToGlobalCCt, mpcEqNums);
     // assemble the global CCt matrix (i.e. add the sub's contributions)
     for(int i=0; i<this->numSubsWithMpcs; ++i)
       this->subsWithMpcs[i]->assembleGlobalCCtsolver(CCtsolver);
   } 
   else { 
     for(int i=0; i<this->numSubsWithMpcs; ++i) 
       this->subsWithMpcs[i]->assembleGlobalCCtsolver(CCtsolver, mpcEqNums);
   }
#ifdef DISTRIBUTED
   CCtsolver->unify(this->fetiCom);
#endif
   //if(this->fetiCom->cpuNum() == 0) ((GenSkyMatrix<Scalar> *)CCtsolver)->print(stderr);
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::factor()
{
  CCtsolver->parallelFactor();
  if(CCtsolver->numRBM() > 0 &&  domain->solInfo().getFetiInfo().contactPrintFlag)
    filePrint(stderr," *** WARNING: Number of singularities in CCt = %d\n", CCtsolver->numRBM());
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::zeroAll()
{
  CCtsolver->zeroAll();
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::reSolve(GenDistrVector<Scalar> &v)
{
  // PJSA: this function is used when applying the generalized preconditioner
  // to the mpc part of the residual
  // computes:  v = (CC^t)^-1 v
  // Step 1. extract mpc part of the residual
  GenVector<Scalar> mpcv1(this->glNumMpc, 0.0);
  execParal2R(this->numSubsWithMpcs, this, &GlobalCCtSolver<Scalar>::extractMpcResidual, v, mpcv1);
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(this->glNumMpc, mpcv1.data());
  for(int i=0; i<this->glNumMpc; ++i)
    if(this->mpcToCpu->num(i) > 1)
      mpcv1[mpcEqNums->firstdof(i)] /= double(this->mpcToCpu->num(i));
#endif
  // Step 2. solve mpcv1 = (CC^t)^-1 mpcv1
  CCtsolver->reSolve(mpcv1);
  // Step 3. insert new mpcv1 back into v
  execParal2R(this->numSubsWithMpcs, this, &GlobalCCtSolver<Scalar>::insertMpcResidual, v, mpcv1);
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::computeSubContributionToGlobalCCt(int i, SimpleNumberer *mpcEqNums)
{
  this->subsWithMpcs[i]->computeSubContributionToGlobalCCt(mpcEqNums);
}

template<class Scalar>
void
GlobalCCtSolver<Scalar>::extractMpcResidual(int iSub, GenDistrVector<Scalar> &v, GenVector<Scalar> &mpcv1)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->extractMpcResidual(subv, mpcv1, mpcEqNums);
}
                                                                                                                                        
template<class Scalar>
void
GlobalCCtSolver<Scalar>::insertMpcResidual(int iSub, GenDistrVector<Scalar> &v, GenVector<Scalar> &mpcv1)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->insertMpcResidual(subv, mpcv1, mpcEqNums);
}


