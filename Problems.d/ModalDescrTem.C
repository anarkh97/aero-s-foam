#include <Problems.d/ModalDescr.h>
#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>
#include <Driver.d/Dynam.h>

template <class Scalar>
ModalDescr<Scalar>::ModalDescr(Domain *d) : ModalBase(d){

  flExchanger = domain->getFileExchanger();

}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::projectForce(Vector &fullF, Vector& modalF){
/*PRE: ModeBase::modesFl have been populated by a call to
         ModeBase::populateFlexModes
 POST: return in modalF, fullF projected onto modesFl
 NOTE: this projection is not valid for displacements, velocities or accelerations
*/
  for(int i = 0; i < numModes; ++i)
    modalF[i] = modesFl[i] * fullF;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::expand(const Vector &modalV, Vector& fullV){
/*PRE: there is at least one modesFl; modesFl are populated
 POST: return in fullV, modalV projected into full space
 NOTE: this expansion is not valid for forces
*/
  fullV.linC(modalV[0], modesFl[0]);
  for(int i = 1; i < numModes; ++i)
    fullV.linAdd(modalV[i], modesFl[i]);
}

//------------------------------------------------------------------------------
template <class Scalar>
void ModalDescr<Scalar>::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
  
}

//------------------------------------------------------------------------------
template <class Scalar>
void ModalDescr<Scalar>::preProcess(){
/*NOTE: call to populateFlexModes gives 1 for 2nd arguement to indicate all
          modes should be read, including those with zero frequency
        see also ModalBase::populateFlexModes
*/
  preProcessBase();
  populateFlexModes(1.0, 1);
  numModes = numFlex;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getTimes(double &dt, double &tmax){

  dt   = domain->solInfo().getTimeStep();
  tmax = domain->solInfo().tmax;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInitState(SysState<Vector> &state){

  initStateBase(state.getDisp(), state.getVeloc(),
    state.getAccel(), state.getPrevVeloc());

/*
  state.getDisp().zero();
  state.getVeloc().zero();
  state.getAccel().zero();
  state.getPrevVeloc().zero();

  // values for state obtained directly from domain
  for(int jj = 0; jj <  domain->numInitDisp(); ++jj)
    state.getDisp()[domain->getInitDisp()[jj].nnum] +=
      domain->getInitDisp()[jj].val;
*/
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getConstForce(Vector &constF){

  domain->computeConstantForce(fullTmpGrav); 
/*
  fullTmpGrav.zero();

  if( domain->gravityFlag()  ) domain->buildGravityForce(fullTmpGrav);
  if( domain->pressureFlag() ) domain->buildPressureForce(fullTmpGrav);
*/
  projectForce(fullTmpGrav, constF);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getSteadyStateParam(int &flag, int &min, int &max, double &tol){

  flag = domain->solInfo().steadyFlag;
  min  = domain->solInfo().steadyMin;
  max  = domain->solInfo().steadyMax;
  tol  = domain->solInfo().steadyTol;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getNewMarkParameters(double &beta, double &gamma,
  double &alphaf, double &alpham){

  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
/*
  fprintf(stderr,"Generalized-Alpha Method: (Newmark)\n");
  fprintf(stderr,"  beta   = %f\n",beta);
  fprintf(stderr,"  gamma  = %f\n",gamma);
  fprintf(stderr,"  alphaf = %f\n",alphaf);
  fprintf(stderr,"  alpham = %f\n",alpham);
*/
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInitialTime(int &tIndex, double &time){

  tIndex = domain->solInfo().initialTimeIndex;
  time   = domain->solInfo().initialTime;
}

//------------------------------------------------------------------------------

template <class Scalar>
double
ModalDescr<Scalar>::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

//------------------------------------------------------------------------------

template <class Scalar>
ModalOps* ModalDescr<Scalar>::buildOps(double mcoef, double ccoef, double kcoef){
//PRE: ModalProbDescr is instantiated; modeData is populated
// POST: operator for time integration
//
  modalOps.M       = new DiagonalMatrix(numModes);
  modalOps.Msolver = new DiagonalMatrix(numModes);
  // see below for damping matrix
  modalOps.K      = new DiagonalMatrix(numModes);
  modalOps.dynMat = new DiagonalMatrix(numModes);

  modalOps.M->setDiag(1.0);
  modalOps.Msolver->setDiag(1.0); // Inverse of M

  int i;
  for(i = 0; i < numModes; ++i){
    (*modalOps.K)[i]      = freqs[i] * freqs[i];
    (*modalOps.dynMat)[i] = kcoef*(*modalOps.K)[i] + mcoef*(*modalOps.M)[i];
  }

  // damping matrix
  double alpha = domain->solInfo().alphaDamp;
  double beta  = domain->solInfo().betaDamp;

  BCond* damping;
  int numDampedModes = geoSource->getModalDamping(damping);

  if(damping){
    modalOps.C = new DiagonalMatrix(numModes);
    modalOps.C->setDiag(0.0);
    for(i = 0; i < numDampedModes; ++i)
      (*modalOps.C)[damping[i].nnum] += 2*damping[i].val*freqs[damping[i].nnum];
    for(i = 0; i < numModes; ++i)
      (*modalOps.dynMat)[i] += ccoef * (*modalOps.C)[i];
  }
  else if(alpha != 0.0 || beta != 0.0){
    modalOps.C = new DiagonalMatrix(numModes);
    for(i = 0; i < numModes; ++i){
      (*modalOps.C)[i] = alpha*(*modalOps.M)[i] + beta*(*modalOps.K)[i];
      (*modalOps.dynMat)[i] += ccoef * (*modalOps.C)[i];
    }
  }
  else{modalOps.C = 0;}

  modalOps.dynMat->invertDiag();

  return (&modalOps);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getQuasiStaticParameters(double &maxVel, double &delta){
  maxVel = domain->solInfo().qsMaxvel;
  delta  = domain->solInfo().delta;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::computeExtForce2(SysState<Vector>& state, Vector &extF,
  Vector &constF, int tIndex, double time, Vector *aeroF,
  double gamma, double alphaf){
/*PRE:
 POST: return in extF, the modalized external force
*/
  domain->template computeExtForce4<double>(fullTmpF, fullTmpGrav, time, 0);
  if(domain->solInfo().aeroFlag >= 0 && tIndex >= 0) {
    domain->buildAeroelasticForce(*aeroF, *prevFrc, tIndex, time, gamma, alphaf);
    fullTmpF += *aeroF;
  }
  projectForce(fullTmpF, extF);
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::getInternalForce(Vector &d, Vector &f, double t){
/*PRE: d is the value of the modal coordinates
 POST: return the modal internal force in f
*/
  for(int i = 0; i < d.size(); ++i)
    f[i] = d[i]*freqs[i]*freqs[i];
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::dynamOutput(int tIndex, ModalOps &ops, Vector &extF,
  Vector *aeroF, SysState<Vector> &state){

  expand(state.getDisp(), fullDsp);
  expand(state.getVeloc(), fullVel);
  expand(state.getAccel(), fullAcc);
  expand(state.getPrevVeloc(), fullPrevVel);

  DynamMat dumDMat;
  domain->dynamOutput(tIndex, bcx, dumDMat, fullTmpF, fullAeroF,
    fullDsp, fullVel, fullAcc, fullPrevVel, vcx);
//  outputModal(state.getDisp(), extF, tIndex);
  outputModal(state, extF, tIndex);
}

//------------------------------------------------------------------------------

template <class Scalar>
int ModalDescr<Scalar>::aeroPreProcess(Vector& d_n, Vector& v_n,
  Vector& a_n, Vector& v_p){
/*PRE: arguements are in modal space
 POST: call domain aeroPreProcess with the expanded state vectors
*/
  expand(d_n, fullDsp);
  expand(v_n, fullVel);
  expand(a_n, fullAcc);
  expand(v_p, fullPrevVel);

  domain->aeroPreProcess(fullDsp, fullVel, fullAcc, fullPrevVel, bcx, vcx);

  return domain->solInfo().aeroFlag;
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::a5TimeLoopCheck(int& parity, double& t, double dt){
/*NOTE: this function copied from SingleDomainDynamic::a5TimeLoopCheck
*/
  if (domain->solInfo().aeroFlag == 5) {
    if (!parity) t -= dt;
    parity = ( parity ? 0 : 1 );
  }
}

//------------------------------------------------------------------------------

template <class Scalar>
void ModalDescr<Scalar>::a5StatusRevise(int parity, SysState<Vector>& curState,
  SysState<Vector>& bkState){
/*NOTE: this function copied from SingleDomainDynamic::a5StatusRevise
*/
  if (domain->solInfo().aeroFlag == 5) {
    if (parity) { // restore

      *prevFrc = *prevFrcBackup;
      curState.getDisp()      = bkState.getDisp();
      curState.getVeloc()     = bkState.getVeloc();
      curState.getAccel()     = bkState.getAccel();
      curState.getPrevVeloc() = bkState.getPrevVeloc();

    } else {      // backup

      *prevFrcBackup = *prevFrc;
      bkState.getDisp()      = curState.getDisp();
      bkState.getVeloc()     = curState.getVeloc();
      bkState.getAccel()     = curState.getAccel();
      bkState.getPrevVeloc() = curState.getPrevVeloc();

    }
  }
}


