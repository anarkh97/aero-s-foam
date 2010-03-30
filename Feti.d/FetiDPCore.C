template<class Scalar> 
void
GenFetiDPSolver<Scalar>::addRstar_gT(int iGroup, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta)
{
#ifdef DISTRIBUTED
  for(int i=0; i<this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup])
      this->sd[i]->addRstar_gT(u.subData(this->sd[i]->localSubNum()), beta);
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->addRstar_gT(u.subData(this->sd[iSub]->localSubNum()), beta);
  }
#endif
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subtractRstar_g(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta)
{
  this->sd[iSub]->subtractRstar_g(u.subData(this->sd[iSub]->localSubNum()), beta);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeProjectedDisplacement(GenDistrVector<Scalar> &u)
{
  int numGtGsing = (GtGtilda) ? GtGtilda->numRBM() : 0;
  if(numGtGsing > 0) {
    // get null space of GtGtilda
    Scalar *zem = new Scalar[numGtGsing*ngrbms];
    GtGtilda->getNullSpace(zem);
    GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
    // build global RBMs 
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::buildGlobalRBMs, X, cornerToSub);

    GenVector<Scalar> beta(numGtGsing, 0.0);
    // 1. compute beta = Rstar_g^t * u
    execParal2R(nGroups1, this, &GenFetiDPSolver<Scalar>::addRstar_gT, u, beta);
#ifdef DISTRIBUTED
    this->fetiCom->globalSum(numGtGsing, beta.data());
#endif
    // build RtR
    GenFullM<Scalar> RtR(numGtGsing, numGtGsing);
    RtR.zero();
    for(int i=0; i<this->nsub; ++i) this->sd[i]->assembleRtR(RtR); // can't be done in parallel
#ifdef DISTRIBUTED
    this->fetiCom->globalSum(numGtGsing*numGtGsing, RtR.getData());
#endif
    RtR.factor();

    // 3. compute beta = (Rstar_g^t Rstar_g)^(-1) * beta
    RtR.reSolve(beta.getData());
    // 4. compute u = u - Rstar_g * beta
    execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::subtractRstar_g, u, beta);
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::projectAlpha(GenVector<Scalar> &alpha, bool flag)
{
  int numGtGsing = (GtGtilda) ? GtGtilda->numRBM() : 0;
  if(numGtGsing > 0) {
    // get null space of GtGtilda
    Scalar *zem = new Scalar[numGtGsing*ngrbms];
    GtGtilda->getNullSpace(zem);
    GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
    // XtX is not necessary if null space has orthogonal basis
    GenFullM<Scalar> XtX(numGtGsing); 
    XtX = X^X; 
    XtX.factor();
    // flag == true: alpha = (I - X*(X^T*X)^{-1}*X^T)*alpha
    // flag == false: alpha = X*(X^T*X)^{-1}*X^T*alpha
    GenVector<Scalar> beta(numGtGsing);
    beta = X^alpha;
    XtX.reSolve(beta.getData());
    if(flag) alpha -= X*beta;
    else alpha = X*beta;
  }  
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeGtG()
{
  filePrint(stderr, " ... Build G matrix and GtG solver  ...\n");
  int i;
  nStatChDual_gtg = nStatChPrimal_gtg = 0;
  startTimerMemory(this->times.coarse1, this->times.memoryGtG);

  // 0. delete previous G if it exists
  deleteG();

  // 1. make local G = Bbar * R
  // this function also initializes GTilda to G
  if(this->fetiInfo->nullSpace == FetiInfo::grbm) // GRBM
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeG);
  else /*if(KccSolver->numRBM() > 0)*/ { // TRBM
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeTrbmG, kccrbms, KccSolver->numRBM(), KccSolver->neqs());
    paralApply(this->nsub, this->sd, &BaseSub::getNumGroupRBM, ngrbmGr);
  }

  // 2. exchange G's between neighboring subdomains
  FSCommPattern<int> *sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
  for(i=0; i<this->nsub; ++i) this->sd[i]->setMpcNeighbCommSize(sPat, 2);
  sPat->finalize();
  paralApply(this->nsub, this->sd, &BaseSub::sendNeighbGrbmInfo, sPat);
  sPat->exchange();  // neighboring subs number of group RBMs and offset
  paralApply(this->nsub, this->sd, &BaseSub::receiveNeighbGrbmInfo, sPat);
  delete sPat;
  // create bodyRbmPat FSCommPattern object, used to send/receive G matrix
  FSCommPattern<Scalar> *gPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU,
                                       FSCommPattern<Scalar>::CopyOnSend, FSCommPattern<Scalar>::NonSym);
  for(i=0; i<this->nsub; ++i) this->sd[i]->setGCommSize(gPat);
  gPat->finalize();
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendG, gPat);
  gPat->exchange();
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::receiveG, gPat);
  delete gPat;

  // 3. build coarse connectivity and equation numberer
  if(this->mpcToSub) {
    Connectivity *mpcToBody = this->mpcToSub->transcon(subToGroup); // PJSA 6-15-06
    Connectivity *bodyToMpc = mpcToBody->reverse();
    Connectivity *bodyToBody_mpc = bodyToMpc->transcon(mpcToBody);
    coarseConnectGtG = bodyToBody_mpc->modify();
    delete mpcToBody; delete bodyToMpc; delete bodyToBody_mpc;
  }
  else {
    Connectivity *bodyToSub = subToBody->reverse();
    coarseConnectGtG = bodyToSub->transcon(subToBody);
    delete bodyToSub;
  }
  eqNumsGtG = new SimpleNumberer(nGroups);
  for(i = 0; i < nGroups; ++i) eqNumsGtG->setWeight(i, ngrbmGr[i]);
  eqNumsGtG->makeOffset();

  // 4. create, assemble and factorize GtG
  GtG = newSolver(this->fetiInfo->auxCoarseSolver, coarseConnectGtG, eqNumsGtG, this->fetiInfo->grbm_tol, GtGsparse);
  //GtG->setPrintNullity(false);
  execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleGtG, 0);
#ifdef DISTRIBUTED
  GtG->unify(this->fetiCom);
#endif
  //cerr << "GtG = \n"; ((SkyMatrix *) GtGsparse)->print();
  startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  GtG->parallelFactor();
  stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  GtGtilda = GtG;

  // 5. check for singularities in GtGstar (representing global RBMs)
  //if(GtG->numRBM() > 0) 
  //  filePrint(stderr, " ... GtG has %d singularities for tol %e ...\n", GtG->numRBM(), this->fetiInfo->grbm_tol);

  stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::makeE(GenDistrVector<Scalar> &f)
{
  // Compute e = R^t * f
  GenVector<Scalar> &e = this->wksp->ret_e();
  e.zero();

  if(this->fetiInfo->nullSpace == FetiInfo::trbm) { // TRBM
    GenVector<Scalar> &fc = this->wksp->ret_fc();
    GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
    Vector dummy(KccSolver->neqs(), 0.0);
    for(int i=0; i<this->nsub; ++i) this->sd[i]->assembleTrbmE(kccrbms, KccSolver->numRBM(), KccSolver->neqs(), e.data(), fr.subData(this->sd[i]->localSubNum()));
#ifdef DISTRIBUTED
    this->fetiCom->globalSum(ngrbms, e.data());
#endif
    for(int i=0; i<KccSolver->numRBM(); ++i)
      for(int j=0; j<KccSolver->neqs(); ++j)
        e[i] += fc[j]*kccrbms[i*KccSolver->neqs()+j];
  } else { // GRBM
    execParal2R(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleE, e, f);
#ifdef DISTRIBUTED
    this->fetiCom->globalSum(ngrbms, e.data());
#endif
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::assembleE(int iGroup, GenVector<Scalar> &e, GenDistrVector<Scalar> &f)
{
#ifdef DISTRIBUTED
  for(int i=0; i<this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup])
      this->sd[i]->assembleE(e, f.subData(this->sd[i]->localSubNum()));
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->assembleE(e, f.subData(this->sd[iSub]->localSubNum()));
  }
#endif
}

template<class Scalar> 
void 
GenFetiDPSolver<Scalar>::addRalpha(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha)
{
  this->sd[iSub]->addRalpha(u.subData(this->sd[iSub]->localSubNum()), alpha);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getRBMs(Scalar *globRBM)
{
  if(GtGtilda) {
    int nRBM = numRBM();
    int iRBM;
    for(iRBM = 0; iRBM < nRBM; ++iRBM) {
      GenStackDistVector<Scalar> R(this->internalDI, globRBM+iRBM*(this->internalDI.len));
      R.zero();
      execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::getGlobalRBM, iRBM, (GenDistrVector<Scalar> &)(R));
    }
  }
  else if(KccSolver) {
    int nc = KccSolver->neqs();
    int nr = numRBM();
    Scalar *R = new Scalar[nr*nc];
    if(nc > 0) KccSolver->getNullSpace(R);
    GenDistrVector<Scalar> vr(internalR);
    int iRBM;
    for(iRBM=0; iRBM<nr; ++iRBM) {
      // note this code should never be required since singularities are always eliminated from Kcc now
      GenStackDistVector<Scalar> v(this->internalDI, globRBM+iRBM*this->internalDI.len);
      GenStackVector<Scalar> vc(nc, R+iRBM*nc);
      vr.zero();
      execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, vr, (GenVector<Scalar> &)(vc));
      execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::KrrReSolve, vr);
      execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, vr, (GenVector<Scalar> &)(vc), 
                  (GenDistrVector<Scalar> &)(v), (GenDistrVector<Scalar> &)(v));  // last argument is a dummy
    }
    delete [] R;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getRBMs(GenDistrVectorSet<Scalar> &globRBM)
{
  if(GtGtilda) {

    int numGtGsing = GtGtilda->numRBM();
    if(numGtGsing > 0 && domain->probType() == SolverInfo::Modal) {
      // get null space of GtGtilda
      Scalar *zem = new Scalar[numGtGsing*ngrbms];
      GtGtilda->getNullSpace(zem);
      GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
      // build global RBMs
      paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::buildGlobalRBMs, X, cornerToSub);
    }

    int nRBM = numRBM();
    int iRBM;
    for(iRBM = 0; iRBM < nRBM; ++iRBM) {
      execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::getGlobalRBM, iRBM, globRBM[iRBM]);
    }
  }
  else if(KccSolver) {
    int nc = KccSolver->neqs();
    int nr = numRBM();
    Scalar *R = new Scalar[nr*nc];
    if(nc > 0) KccSolver->getNullSpace(R);
    GenDistrVector<Scalar> vr(internalR);
    int iRBM;
    for(iRBM=0; iRBM<nr; ++iRBM) {
      // note this code should never be required since singularities are always eliminated from Kcc now
      GenStackVector<Scalar> vc(nc, R+iRBM*nc);
      vr.zero();
      execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, vr, (GenVector<Scalar> &)(vc));
      execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::KrrReSolve, vr);
      execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, vr, (GenVector<Scalar> &)(vc),
                  globRBM[iRBM],globRBM[iRBM]); // last argument is a dummy
    }
    if(R) delete [] R;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R)
{
  Scalar *localRvec = R.subData(this->sd[iSub]->localSubNum());
  this->sd[iSub]->getGlobalRBM(iRBM, localRvec);
}

template<>
inline void
GenFetiDPSolver<DComplex>::split(int iSub, GenDistrVector<DComplex> &v, GenDistrVector<DComplex> &v_f, GenDistrVector<DComplex> &v_c,
                                 GenDistrVector<DComplex> &v_p)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::split(...) is not implemented \n");
}

template<>
inline void
GenFetiDPSolver<double>::split(int iSub, GenDistrVector<double> &v, GenDistrVector<double> &v_f, GenDistrVector<double> &v_c,
                               GenDistrVector<double> &v_p)
{
  this->sd[iSub]->split(v.subData(this->sd[iSub]->localSubNum()), v_f.subData(this->sd[iSub]->localSubNum()), v_c.subData(this->sd[iSub]->localSubNum()),
                  v_p.subData(this->sd[iSub]->localSubNum()));
}

template<>
inline void
GenFetiDPSolver<DComplex>::chop(int iSub, GenDistrVector<DComplex> &v, GenDistrVector<DComplex> &v_c, double tol, int chop_flag)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::chop(...) is not implemented \n");
}

template<>
inline void
GenFetiDPSolver<double>::chop(int iSub, GenDistrVector<double> &v, GenDistrVector<double> &v_c, double tol, int chop_flag)
{
  this->sd[iSub]->chop(v.subData(this->sd[iSub]->localSubNum()), v_c.subData(this->sd[iSub]->localSubNum()), tol, chop_flag);
}
