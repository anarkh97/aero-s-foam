#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include "DistrGalerkinProjectionSolver.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"

#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <limits>

namespace Rom {

DistrExplicitPodPostProcessor::DistrExplicitPodPostProcessor(DecDomain *d, StaticTimers* _times, DistrGeomState *_geomState = 0, Corotator ***_allCorot = 0) :
    MultiDomDynPostProcessor(d, _times, _geomState, _allCorot)
{

      decDomain = d;
      geomState = _geomState;
      times = _times;
      oinfo = geoSource->getOutputInfo();

      numOutInfo = geoSource->getNumOutInfo();


      for (int iOut = 0; iOut < numOutInfo; iOut++) {
       switch(oinfo[iOut].type){
         case OutputInfo::Accel6 : 
           oinfo[iOut].filptr = fopen(oinfo[iOut].filename,"wb");
           filePrint( oinfo[iOut].filptr, "0\n"); 
           break;
         case OutputInfo::Disp6DOF :
           oinfo[iOut].filptr = fopen(oinfo[iOut].filename,"wb");
           filePrint( oinfo[iOut].filptr, "1\n");
           break;
         case OutputInfo::Velocity6 : 
           oinfo[iOut].filptr = fopen(oinfo[iOut].filename,"wb");
           filePrint( oinfo[iOut].filptr, "2\n"); 
           break; 
         default:
           filePrint(stderr, " ...ROM output only supports Acceleration, Displacement, and Velocity... \n");
           filePrint(stderr, " output type selected is %d \n", oinfo[iOut].type);
       }
      }

}

DistrExplicitPodPostProcessor::~DistrExplicitPodPostProcessor() {

/*  for (int iOut = 0; iOut < numOutInfo; iOut++) {
           if(oinfo[iOut].filptr) fclose(oinfo[iOut].filptr);
      }
*/
}

void
DistrExplicitPodPostProcessor::printPODSize(int PODsize) {

    podSize = PODsize;

    for (int iOut = 0; iOut < numOutInfo; iOut++) {
     switch(oinfo[iOut].type){
       case OutputInfo::Accel6 :
         filePrint( oinfo[iOut].filptr, "%d\n", PODsize);
         break;
       case OutputInfo::Disp6DOF :
         filePrint( oinfo[iOut].filptr, "%d\n", PODsize);
         break;
       case OutputInfo::Velocity6 :
         filePrint( oinfo[iOut].filptr, "%d\n", PODsize);
         break;
       default:
         break;
     }
    }

}

void
DistrExplicitPodPostProcessor::dynamOutput(int tIndex, double t, MDDynamMat &dynOps, DistrVector &distForce,
                                                          DistrVector *distAeroF, SysState<DistrVector>& distState) {


  //all MPI processes have a full copy of reduced coordinates, only master processes needs to print
  //valgrind shows uninitialized conditionals junps here, need to figure out why
  int p = std::numeric_limits<double>::digits10+1;
  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(tIndex % oinfo[iOut].interval == 0) {

      switch(oinfo[iOut].type){
         case OutputInfo::Accel6 :
           {
           filePrint(oinfo[iOut].filptr, "%.*e\n", p, t); // print timestamp
           for(int i = 0; i<podSize; i++) {
             filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getAccel()[i]);
           }
           filePrint(oinfo[iOut].filptr, "\n");}
           break;
         case OutputInfo::Disp6DOF :
           {
           filePrint(oinfo[iOut].filptr, "%.*e\n", p, t); // print timestamp
           for(int i = 0; i<podSize; i++) {
             filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getDisp()[i]);
           }
           filePrint(oinfo[iOut].filptr, "\n");}
           break;
         case OutputInfo::Velocity6 : 
           {
           filePrint(oinfo[iOut].filptr, "%.*e\n", p, t); // print timestamp
           for(int i = 0; i<podSize; i++) {
             filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getVeloc()[i]);
           }
           filePrint(oinfo[iOut].filptr, "\n");}
           break;
         default:
           filePrint(stderr, " ...ROM output only supports Acceraltion, Displacement, and Velocity... \n");
      }
    }
  }
}

DistrExplicitPodPostProcessor *
DistrExplicitPodProjectionNonLinDynamicBase::getPostProcessor() {

   mddPostPro = new DistrExplicitPodPostProcessor(decDomain, times, geomState, allCorot);
   mddPostPro->printPODSize(projectionBasis_.numVectors());

   return mddPostPro;

}

DistrExplicitPodProjectionNonLinDynamicBase::DistrExplicitPodProjectionNonLinDynamicBase(Domain *_domain) :
  MultiDomainDynam(_domain),
  domain(_domain)
{}

void
DistrExplicitPodProjectionNonLinDynamicBase::preProcess() {
  {MultiDomainDynam::preProcess();
  //preProcessing for reduced order basis/////////////////////////////////////////////////
  FileNameInfo fileInfo; 
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  projectionBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = projectionBasis_.begin(),
                               it_end = projectionBasis_.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);
    
    podBasisFile.currentStateIndexInc();
  }}

  ///////////////////////////////////////////////////////////////////////////////////////

  //preProcessing for solution vecotor information///////////////////////////////////////
  //each subdomain gets full copy of reduced coordinates
  {reducedInfo.domLen = new int[MultiDomainDynam::solVecInfo().numDom]; 
  reducedInfo.numDom = MultiDomainDynam::solVecInfo().numDom;
  int totLen = 0;
  for(int iSub = 0; iSub < MultiDomainDynam::solVecInfo().numDom; ++iSub) {
    reducedInfo.domLen[iSub] = projectionBasis_.numVec();
    totLen += reducedInfo.domLen[iSub];
  }

  reducedInfo.len = totLen;
  reducedInfo.setMasterFlag(); //not correct masterflag, but doesn't matter unless computing norm
  }
  //////////////////////////////////////////////////////////////////////////////////////
  
  //preProcessing for dummy working variables///////////////////////////////////////////////////
  {fExt      = new DistrVector(MultiDomainDynam::solVecInfo());
  fInt       = new DistrVector(MultiDomainDynam::solVecInfo());
  cnst_fBig  = new DistrVector(MultiDomainDynam::solVecInfo());
  aero_fBig  = new DistrVector(MultiDomainDynam::solVecInfo());
  d_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  v_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  a_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  v_p        = new DistrVector(MultiDomainDynam::solVecInfo());
  tempVec    = new DistrVector(MultiDomainDynam::solVecInfo());
  dummyState = new SysState<DistrVector>(*d_n, *v_n, *a_n, *v_p);
  times      = new StaticTimers;}
  ///////////////////////////////////////////////////////////////////////////////////////

}

const DistrInfo &
DistrExplicitPodProjectionNonLinDynamicBase::solVecInfo() {
  return reducedVecInfo();
}

DistrInfo &
DistrExplicitPodProjectionNonLinDynamicBase::reducedVecInfo() {
  return reducedInfo; //ROB size
}

void 
DistrExplicitPodProjectionNonLinDynamicBase::printFullNorm(DistrVector &v) {

  normalizedBasis_.projectUp(v,*tempVec);

  filePrint(stderr,"%1.4e\n",tempVec->norm());

}

void
DistrExplicitPodProjectionNonLinDynamicBase::getInitState(SysState<DistrVector> & _curState) {
  //project initial state into reduced coordinates

  MultiDomainDynam::getInitState( *dummyState );

  DistrVector &_d_n = _curState.getDisp(); 
  DistrVector &_v_n = _curState.getVeloc();
  DistrVector &_a_n = _curState.getAccel();
  DistrVector &_v_p = _curState.getPrevVeloc();

  //this projection doesn't do anything since the projectionBasis_ isn't initialized yet
  //this only matters if we have a case where the initial conditions are other than 0
  //need to fixe this if we want to use resart
  projectionBasis_.projectDown( *d_n, _d_n);
  projectionBasis_.projectDown( *v_n, _v_n);
  projectionBasis_.projectDown( *a_n, _a_n);
  projectionBasis_.projectDown( *v_p, _v_p);
}

void 
DistrExplicitPodProjectionNonLinDynamicBase::updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n1) {
  //update geomState for Fint, but no need to update displacment vector from geometry 

  DistrVector temp1(solVecInfo());
  temp1 = dt_n_h*v_n_h;

  normalizedBasis_.projectUp( temp1, *d_n); 

  geomState->update(*d_n, 2);

  d_n1 += temp1;  //we save the increment vectors for postprocessing

  bool updateVelocityInGeomState = true; // TODO only set this to true when necessary
  if(updateVelocityInGeomState) {
    normalizedBasis_.projectUp(v_n_h, *v_n);
    geomState->setVelocity(*v_n, 2);
  }

}

void DistrExplicitPodProjectionNonLinDynamicBase::getConstForce(DistrVector& v)
{
  //we really don't need to project down here since cnst_fBig is stored inside the probDesc class
  //just a formality. 
  MultiDomainDynam::getConstForce(*cnst_fBig);
  normalizedBasis_.projectDown(*cnst_fBig,v);
  cnst_fBig->zero();
}

void
DistrExplicitPodProjectionNonLinDynamicBase::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {
  //Build internal force and project into reduced coordinates

  MultiDomainDynam::getInternalForce( *d_n, *fInt, t, tIndex);
  //compute residual here to prevent having to project into reduced basis twice
  *a_n = *fInt - *fExt;

  bool hasRot = true; // TODO: only do this when model has rotation dofs see also buildOps
  if(hasRot) {
    geomState->transform(*a_n, 3);
    fullMassSolver->reSolve(*a_n);
    geomState->transform(*a_n, 2);
    DistrVector toto(*a_n);
    dynMat->M->mult(toto, *a_n);
  }

  normalizedBasis_.projectDown(*a_n,f); 
}

void
DistrExplicitPodProjectionNonLinDynamicBase::computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f, int tIndex,
                        double t, DistrVector *aero_f,
                        double gamma, double alphaf) {

  f = cnst_f;
  MultiDomainDynam::computeExtForce2( *dummyState, *fExt, *cnst_fBig, tIndex, t, aero_fBig, gamma, alphaf);

  //f += cnst_f; should implement another version were the constant force is added 
  // in the reduced coordinates to get a little more speed
}

MDDynamMat *
DistrExplicitPodProjectionNonLinDynamicBase::buildOps(double mCoef, double cCoef, double kCoef) {
  MDDynamMat *result = MultiDomainDynam::buildOps(mCoef, cCoef, kCoef);
  assert(result->M);

  const GenSubDOp<double> &fullMass = *(result->M);
  renormalized_basis(fullMass, projectionBasis_, normalizedBasis_);

  std::auto_ptr<DistrGalerkinProjectionSolver> solver(new DistrGalerkinProjectionSolver(normalizedBasis_));

  bool hasRot = true; // TODO: only do this when model has rotation dofs see also getInternalForce
  if(hasRot) {
    fullMassSolver = result->dynMat;
  }
  else {
    delete result->dynMat;
  }

  result->dynMat = solver.release();

  return result;
}

} // end namespace Rom
