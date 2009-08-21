#include <Feti.d/Feti.h>

template<> 
void
GenFetiDPSolver<double>::addRstar_gT(int iGroup, GenDistrVector<double> &u, GenVector<double> &beta)
{
#ifdef DISTRIBUTED
  for(int i=0; i<nsub; ++i) {
    if(sd[i]->group == groups[iGroup])
      sd[i]->addRstar_gT(u.subData(sd[i]->localSubNum()), beta);
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    sd[iSub]->addRstar_gT(u.subData(sd[iSub]->localSubNum()), beta);
  }
#endif
}

template<> 
void
GenFetiDPSolver<double>::subtractRstar_g(int iSub, GenDistrVector<double> &u, GenVector<double> &beta)
{
  sd[iSub]->subtractRstar_g(u.subData(sd[iSub]->localSubNum()), beta);
}

template<>
void
GenFetiDPSolver<DComplex>::computeProjectedDisplacement(GenDistrVector<DComplex> &u)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::computeProjectedDisplacement(...) not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::computeProjectedDisplacement(GenDistrVector<double> &u)
{
  int numGtGsing = (GtGtilda) ? GtGtilda->numRBM() : 0;
  if(numGtGsing > 0) {
    // get null space of GtGtilda
    double *zem = new double[numGtGsing*ngrbms];
    GtGtilda->getNullSpace(zem);
    FullM X(zem, ngrbms, numGtGsing, 1);
    // build global RBMs 
    paralApply(nsub, sd, &BaseSub::buildGlobalRBMs, X, cornerToSub);

    GenVector<double> beta(numGtGsing, 0.0);
    // 1. compute beta = Rstar_g^t * u
    execParal2R(nGroups1, this, &GenFetiDPSolver<double>::addRstar_gT, u, beta);
#ifdef DISTRIBUTED
    fetiCom->globalSum(numGtGsing, beta.data());
#endif
    // build RtR
    FullM RtR(numGtGsing, numGtGsing);
    RtR.zero();
    for(int i=0; i<nsub; ++i) sd[i]->assembleRtR(RtR);  // can't be done in parallel
#ifdef DISTRIBUTED
    fetiCom->globalSum(numGtGsing*numGtGsing, RtR.getData());
#endif
    RtR.factor();

    // 3. compute beta = (Rstar_g^t Rstar_g)^(-1) * beta
    RtR.reSolve(beta.getData());
    // 4. compute u = u - Rstar_g * beta
    execParal2R(nsub, this, &GenFetiDPSolver<double>::subtractRstar_g, u, beta);
  }
}

template<>
bool
GenFetiDPSolver<DComplex>::inconsistent(GenVector<DComplex> &e, bool print_warning)
{
  if(ngrbms) filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::inconsistent(...) not implemented \n");
  return false;
}

template<>
bool
GenFetiDPSolver<double>::inconsistent(GenVector<double> &e, bool print_warning)
{
  // check X^t * e = 0, where X is the null space of Gtilda^T*Gtilda
  if(GtGtilda->numRBM() == 0) return false;

  double *zem = new double[GtGtilda->numRBM()*ngrbms];
  GtGtilda->getNullSpace(zem);
  FullM X(zem, ngrbms, GtGtilda->numRBM(), 1);

  GenVector<double> q = X ^ e;
  double inc_err = 0.0;
  bool inc = ((inc_err = q.norm()) >= fetiInfo->equi_tol);
  if(inc && print_warning && myCPU == 0) cerr << "warning: inconsistent constraints in dual-active set (error = " << inc_err << ", equi_tol = " << fetiInfo->equi_tol << ")\n";
  return inc;
}

template<>
void
GenFetiDPSolver<DComplex>::makeGtG()
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::makeGtG() is not implemented");
}

template<> 
void
GenFetiDPSolver<double>::makeGtG()
{
  filePrint(stderr, " ... Build G matrix and GtG solver  ...\n");
  int i;
  nStatChDual_gtg = nStatChPrimal_gtg = 0;
  startTimerMemory(times.coarse1, times.memoryGtG);

  // 1. make local G = Bbar * R
  // this function also initializes GTilda to G
  paralApply(nsub, sd, &GenSubDomain<double>::makeG);

  // 2. exchange G's between neighboring subdomains
  FSCommPattern<int> *sPat = new FSCommPattern<int>(fetiCom, cpuToSub, myCPU, FSCommPattern<int>::CopyOnSend);
  for(i=0; i<nsub; ++i) sd[i]->setMpcNeighbCommSize(sPat, 2);
  sPat->finalize();
  paralApply(nsub, sd, &BaseSub::sendNeighbGrbmInfo, sPat);
  sPat->exchange();  // neighboring subs number of group RBMs and offset
  paralApply(nsub, sd, &BaseSub::receiveNeighbGrbmInfo, sPat);
  delete sPat;
  // create bodyRbmPat FSCommPattern object, used to send/receive rigid body modes
  FSCommPattern<double> *grbmPat = new FSCommPattern<double>(fetiCom, cpuToSub, myCPU,
                                       FSCommPattern<double>::CopyOnSend, FSCommPattern<double>::NonSym);
  for(i=0; i<nsub; ++i) sd[i]->setGrbmCommSize(grbmPat);
  grbmPat->finalize();
  paralApply(nsub, sd, &BaseSub::sendG, grbmPat);
  grbmPat->exchange();
  paralApply(nsub, sd, &BaseSub::receiveG, grbmPat);
  delete grbmPat;

  // 3. build coarse connectivity and equation numberer
  if(mpcToSub) {
    Connectivity *mpcToBody = mpcToSub->transcon(subToGroup); // PJSA 6-15-06
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
  if(fetiInfo->outerloop == FetiInfo::CGAL && fetiInfo->auxCoarseSolver != FetiInfo::skyline) {
    fetiInfo->auxCoarseSolver = FetiInfo::skyline;
    if(myCPU == 0) cerr << " *** WARNING: selected aux_coarse_sover type not supported for outerloop CGAL, switching to skyline\n";
  }
  GtG = newSolver(fetiInfo->auxCoarseSolver, coarseConnectGtG, eqNumsGtG, fetiInfo->grbm_tol, GtGsparse);
  GtG->setPrintNullity(false);
  execParal(nGroups1, this, &GenFetiDPSolver<double>::assembleGtG, 0);
#ifdef DISTRIBUTED
  GtG->unify(fetiCom);
#endif
  startTimerMemory(times.pfactor, times.memoryGtGsky);
  GtG->parallelFactor();
  stopTimerMemory(times.pfactor, times.memoryGtGsky);
  GtGtilda = GtG;

  // 5. check for singularities in GtGstar (representing global RBMs)
  if(GtG->numRBM() > 0) 
    filePrint(stderr, " ... GtG has %d singularities for tol %e ...\n", GtG->numRBM(), fetiInfo->grbm_tol);

  stopTimerMemory(times.coarse1, times.memoryGtG);
}

template<>
void
GenFetiDPSolver<DComplex>::makeE(GenDistrVector<DComplex> &f)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::makeE(...) not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::makeE(GenDistrVector<double> &f)
{
  // Compute e = R^t * f
  GenVector<double> &e = this->wksp->ret_e();
  e.zero();
  execParal2R(nGroups1, this, &GenFetiDPSolver<double>::assembleE, e, f);
#ifdef DISTRIBUTED
  fetiCom->globalSum(ngrbms, e.data());
#endif
  e_copy = e; // keep a copy of the original rhs, used in project(...)
  if(this->fetiInfo->outerloop == FetiInfo::CGAL) GtG->forward(e); 
  ee = e*e;
  if(ee != 0.0) this->fetiInfo->equi_tol *= sqrt(ee); // convert equi_tol to absolute tolerance
  else ee = 1.0;
}

template<>
void
GenFetiDPSolver<double>::assembleE(int iGroup, GenVector<double> &e, GenDistrVector<double> &f)
{
#ifdef DISTRIBUTED
  for(int i=0; i<nsub; ++i) {
    if(sd[i]->group == groups[iGroup])
      sd[i]->assembleE(e, f.subData(sd[i]->localSubNum()));
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    sd[iSub]->assembleE(e, f.subData(sd[iSub]->localSubNum()));
  }
#endif
}

template<>
void
GenFetiDPSolver<DComplex>::addRalpha(int iSub, GenDistrVector<DComplex> &u, GenVector<DComplex> &alpha)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::addRalpha(...) not implemented \n");
}

template<> 
void 
GenFetiDPSolver<double>::addRalpha(int iSub, GenDistrVector<double> &u, GenVector<double> &alpha)
{
  sd[iSub]->addRalpha(u.subData(sd[iSub]->localSubNum()), alpha);
}

template<>
void
GenFetiDPSolver<DComplex>::getRBMs(DComplex *globRBM)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::getRBMs(DComplex *globRBM) is not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::getRBMs(double *globRBM)
{
  if(GtGtilda) {
    int nRBM = numRBM();
    int iRBM;
    for(iRBM = 0; iRBM < nRBM; ++iRBM) {
      GenStackDistVector<double> R(internalDI, globRBM+iRBM*(internalDI.len));
      R.zero();
      execParal2R(nsub, this, &GenFetiDPSolver<double>::getGlobalRBM, iRBM, (GenDistrVector<double> &)(R));
    }
  }
  else if(KccSolver) {
    int nc = KccSolver->neqs();
    int nr = numRBM();
    double *R = new double[nr*nc];
    if(nc > 0) KccSolver->getRBMs(R);
    GenDistrVector<double> vr(internalR);
    int iRBM;
    for(iRBM=0; iRBM<nr; ++iRBM) {
      // note this code should never be required since singularities are always eliminated from Kcc now
      GenStackDistVector<double> v(internalDI, globRBM+iRBM*internalDI.len);
      GenStackVector<double> vc(nc, R+iRBM*nc);
      vr.zero();
      execParal2R(nsub, this, &GenFetiDPSolver<double>::multKrc, vr, (GenVector<double> &)(vc));
      execParal1R(nsub, this, &GenFetiDPSolver<double>::KrrReSolve, vr);
      execParal4R(nsub, this, &GenFetiDPSolver<double>::mergeUr, vr, (GenVector<double> &)(vc), 
                  (GenDistrVector<double> &)(v), (GenDistrVector<double> &)(v));  // last argument is a dummy
    }
    delete [] R;
  }
}

template<>
void
GenFetiDPSolver<DComplex>::getRBMs(GenDistrVectorSet<DComplex> &globRBM)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::getRBMs(...) not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::getRBMs(GenDistrVectorSet<double> &globRBM)
{
  if(GtGtilda) {

    int numGtGsing = GtGtilda->numRBM();
    if(numGtGsing > 0 && domain->probType() == SolverInfo::Modal) {
      // get null space of GtGtilda
      double *zem = new double[numGtGsing*ngrbms];
      GtGtilda->getNullSpace(zem);
      FullM X(zem, ngrbms, numGtGsing, 1);
      // build global RBMs
      paralApply(nsub, sd, &BaseSub::buildGlobalRBMs, X, cornerToSub);
    }

    int nRBM = numRBM();
    int iRBM;
    for(iRBM = 0; iRBM < nRBM; ++iRBM) {
      execParal2R(nsub, this, &GenFetiDPSolver<double>::getGlobalRBM, iRBM, globRBM[iRBM]);
    }
  }
  else if(KccSolver) {
    int nc = KccSolver->neqs();
    int nr = numRBM();
    double *R = new double[nr*nc];
    if(nc > 0) KccSolver->getRBMs(R);
    GenDistrVector<double> vr(internalR);
    int iRBM;
    for(iRBM=0; iRBM<nr; ++iRBM) {
      // note this code should never be required since singularities are always eliminated from Kcc now
      GenStackVector<double> vc(nc, R+iRBM*nc);
      vr.zero();
      execParal2R(nsub, this, &GenFetiDPSolver<double>::multKrc, vr, (GenVector<double> &)(vc));
      execParal1R(nsub, this, &GenFetiDPSolver<double>::KrrReSolve, vr);
      execParal4R(nsub, this, &GenFetiDPSolver<double>::mergeUr, vr, (GenVector<double> &)(vc),
                  globRBM[iRBM],globRBM[iRBM]); // last argument is a dummy
    }
    if(R) delete [] R;
  }
}

template<>
void
GenFetiDPSolver<DComplex>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<DComplex> &R)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::getGlobalRBM(...) not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<double> &R)
{
  double *localRvec = R.subData(sd[iSub]->localSubNum());
  sd[iSub]->getGlobalRBM(iRBM, localRvec);
}

template<>
void
GenFetiDPSolver<DComplex>::split(int iSub, GenDistrVector<DComplex> &v, GenDistrVector<DComplex> &v_f, GenDistrVector<DComplex> &v_c,
                                 GenDistrVector<DComplex> &v_p)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::split(...) is not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::split(int iSub, GenDistrVector<double> &v, GenDistrVector<double> &v_f, GenDistrVector<double> &v_c,
                               GenDistrVector<double> &v_p)
{
  sd[iSub]->split(v.subData(sd[iSub]->localSubNum()), v_f.subData(sd[iSub]->localSubNum()), v_c.subData(sd[iSub]->localSubNum()),
                  v_p.subData(sd[iSub]->localSubNum()));
}

template<>
void
GenFetiDPSolver<DComplex>::chop(int iSub, GenDistrVector<DComplex> &v, GenDistrVector<DComplex> &v_c, double tol, int chop_flag)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::chop(...) is not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::chop(int iSub, GenDistrVector<double> &v, GenDistrVector<double> &v_c, double tol, int chop_flag)
{
  sd[iSub]->chop(v.subData(sd[iSub]->localSubNum()), v_c.subData(sd[iSub]->localSubNum()), tol, chop_flag);
}

template<>
void
GenFetiDPSolver<DComplex>::subQuotient(int iSub, GenDistrVector<DComplex> &q, GenDistrVector<DComplex> &lambda, GenDistrVector<DComplex> &p)
{
  filePrint(stderr, " *** WARNING: GenFetiDPSolver<DComplex>::subQuotient(...) is not implemented \n");
}

template<>
void
GenFetiDPSolver<double>::subQuotient(int iSub, GenDistrVector<double> &q, GenDistrVector<double> &lambda, GenDistrVector<double> &p)
{
  sd[iSub]->quotient(q.subData(sd[iSub]->localSubNum()), lambda.subData(sd[iSub]->localSubNum()), p.subData(sd[iSub]->localSubNum()));
}

