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

#include <Threads.d/PHelper.h>

#include <Driver.d/DecDomain.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <limits>

namespace Rom {

DistrExplicitPodPostProcessor::DistrExplicitPodPostProcessor(DecDomain *d, StaticTimers* _times, DistrGeomState *_geomState = 0, Corotator ***_allCorot = 0) :
    MultiDomDynPostProcessor(d, _times, _geomState, _allCorot),
    DispSensorValues(),
    AccSensorValues(),
    VelSensorValues(),
    all_cdsa(NULL)
{
  decDomain = d;
  geomState = _geomState;
  times = _times;
  oinfo = geoSource->getOutputInfo();

  numOutInfo = geoSource->getNumOutInfo();

  DispSensor = false;
  AccSensor  = false;
  VelSensor  = false;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber != -1) {
          AccSensor = true;
        } else {
          if(!structCom || structCom->myID() == 0) oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "0\n"); 
        }
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber != -1) {
          DispSensor = true;
        } else {
          if(!structCom || structCom->myID() == 0) oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "1\n");
        }
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber != -1) {
          VelSensor = true;
        } else {
          if(!structCom || structCom->myID() == 0) oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          filePrint(oinfo[iOut].filptr, "2\n");
        }
        break; 
      default:
        filePrint(stderr, " ... ROM output only supports Acceleration, Displacement, and Velocity ...\n");
        filePrint(stderr, "     output type selected is %d \n", oinfo[iOut].type);
    }
  }

  if(DispSensor || VelSensor || AccSensor) {
    for(int iSub = 0; iSub < decDomain->getNumSub(); iSub++) {
      int nOutNodes = decDomain->getSubDomain(iSub)->getNumNodalOutput();
      if(nOutNodes) {
        geoSource->openOutputFiles(decDomain->getSubDomain(iSub)->getOutputNodes(),
                                   decDomain->getSubDomain(iSub)->getOutIndex(), nOutNodes);
      }
    }
  }

  if(structCom) structCom->sync();
  nodeVector.resize(decDomain->getNumSub());
  execParal(decDomain->getNumSub(), this, &DistrExplicitPodPostProcessor::subBuildSensorNodeVector);
}

DistrExplicitPodPostProcessor::~DistrExplicitPodPostProcessor() {

  if(all_cdsa) delete [] all_cdsa;
}

void
DistrExplicitPodPostProcessor::printPODSize(int PODsize) {

  podSize = PODsize;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      default:
        break;
    }
  }
}

void
DistrExplicitPodPostProcessor::makeSensorBasis(DistrVecBasis *fullBasis) {

  all_cdsa = new DofSetArray * [decDomain->getNumSub()];

  for(int i = 0; i < decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();

  SensorBasis.dimensionIs(fullBasis->numVec(), fullBasis->vectorInfo());
  SensorBasis = *fullBasis;
  SensorBasis.makeSparseBasis(nodeVector,all_cdsa);

  // allocate space for containers to hold sensor values
  if(DispSensor)
    new (&DispSensorValues) GenDistrVector<double>(SensorBasis.vectorInfo());
  if(AccSensor)
    new (&AccSensorValues) GenDistrVector<double>(SensorBasis.vectorInfo());
  if(VelSensor)
    new (&VelSensorValues) GenDistrVector<double>(SensorBasis.vectorInfo());
}
  
void
DistrExplicitPodPostProcessor::subBuildSensorNodeVector(int iSub) {

 std::vector<int> &subSensorNodes = nodeVector[iSub];

 //load vector of sensor nodes converted to local numbering
 for (int iOut = 0; iOut < numOutInfo; iOut++) {
   if(oinfo[iOut].nodeNumber != -1) {
     int locNode = decDomain->getSubDomain(iSub)->globalToLocal(oinfo[iOut].nodeNumber);
     if(locNode > -1) 
       subSensorNodes.push_back(locNode);
   }
 }

 //if multiple outputs are requested for a single node, cull redundancies from node vector
 std::sort(subSensorNodes.begin(), subSensorNodes.end());
 std::vector<int>::iterator packedNodeIt = std::unique(subSensorNodes.begin(), subSensorNodes.end());
 subSensorNodes.resize(packedNodeIt-subSensorNodes.begin());
}

void
DistrExplicitPodPostProcessor::subPrintSensorValues(int iSub, GenDistrVector<double> &SensorData, OutputInfo *OINFO, double *time) {

  // XXX the rotation vector should be renormalized, and the angular velocity/acceleration should be transformed to convected. 
  int locNode = decDomain->getSubDomain(iSub)->globalToLocal(OINFO->nodeNumber);
  if(locNode > -1) { // if node is -1, sensor node is not in this subdomain, don't print
    for(int k = 0; k < decDomain->getSubDomain(iSub)->getNumNodalOutput(); ++k) { // not the best way
      if(locNode == decDomain->getSubDomain(iSub)->getOutputNodes()[k]) {

        int ndofs;
        int dofs[6];
        if(OINFO->type == OutputInfo::Disp6DOF || OINFO->type == OutputInfo::Velocity6 || OINFO->type == OutputInfo::Accel6) {
          ndofs = all_cdsa[iSub]->number(locNode, DofSet::XYZdisp | DofSet::XYZrot, dofs);
        }
        else {
          ndofs = all_cdsa[iSub]->number(locNode, DofSet::XYZdisp, dofs);
        }
        int w = OINFO->width;
        int p = OINFO->precision; 
        fprintf(OINFO->filptr, "  % *.*E  ", w, p, *time);
        for(int j = 0; j < ndofs; ++j) {
          if(dofs[j] != -1) {
            fprintf(OINFO->filptr, " % *.*E", w, p, SensorData.subData(iSub)[dofs[j]]);
          }
          else {
            fprintf(OINFO->filptr, " % *.*E", w, p, 0.0); // TODO constrained dofs
          }
        }
        fprintf(OINFO->filptr, "\n");
        fflush(OINFO->filptr);
        break;
      }
    }
  }
}

void
DistrExplicitPodPostProcessor::dynamOutput(int tIndex, double t, MDDynamMat &dynOps, DistrVector &distForce,
                                           DistrVector *distAeroF, SysState<DistrVector>& distState) {

  //all MPI processes have a full copy of reduced coordinates, only master processes needs to print
  int p = std::numeric_limits<double>::digits10+1;

  bool DispProjected = false, AccProjected = false, VelProjected = false;

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(tIndex % oinfo[iOut].interval == 0) {

      switch(oinfo[iOut].type) {
         case OutputInfo::Accel6 : case OutputInfo::Acceleration :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getAccel()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             } else {
               if(!AccProjected) {
                 SensorBasis.projectUp(distState.getAccel(), AccSensorValues);
                 AccProjected = true;
               }
               execParal3R(decDomain->getNumSub(), this, &DistrExplicitPodPostProcessor::subPrintSensorValues, AccSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         case OutputInfo::Disp6DOF : case OutputInfo::Displacement : 
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getDisp()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             } else {
               if(!DispProjected) {
                 SensorBasis.projectUp(distState.getDisp(), DispSensorValues);
                 DispProjected = true;
               }
               execParal3R(decDomain->getNumSub(), this, &DistrExplicitPodPostProcessor::subPrintSensorValues, DispSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         case OutputInfo::Velocity6 : case OutputInfo::Velocity :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getVeloc()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             } else {
               if(!VelProjected) {
                 SensorBasis.projectUp(distState.getVeloc(), VelSensorValues);
                 VelProjected = true;
               }
               execParal3R(decDomain->getNumSub(), this, &DistrExplicitPodPostProcessor::subPrintSensorValues, VelSensorValues, &oinfo[iOut], &t);
             }
           }
           break;
         default:
           filePrint(stderr, " ... ROM output only supports Acceleration, Displacement, and Velocity ... \n");
      }
    }
  }
}

DistrExplicitPodPostProcessor *
DistrExplicitPodProjectionNonLinDynamicBase::getPostProcessor() {

  mddPostPro = new DistrExplicitPodPostProcessor(decDomain, times, geomState, allCorot);
  mddPostPro->printPODSize(normalizedBasis_.numVectors());
  mddPostPro->makeSensorBasis(&normalizedBasis_);

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
  std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
  fileName.append(".normalized");

  DistrBasisInputFile BasisFile(fileName); //read in mass-normalized basis
  filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());
  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, BasisFile.stateCount()) :
                                     BasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  normalizedBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo()); 

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = normalizedBasis_.begin(),
                               it_end = normalizedBasis_.end();
                               it != it_end; ++it) {
    assert(BasisFile.validCurrentState());

    BasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);
    
    BasisFile.currentStateIndexInc();
  }}

  ///////////////////////////////////////////////////////////////////////////////////////

  //preProcessing for solution vecotor information///////////////////////////////////////
  {reducedInfo.domLen = new int[MultiDomainDynam::solVecInfo().numDom]; 
  reducedInfo.numDom = MultiDomainDynam::solVecInfo().numDom;
  int totLen = 0;
  for(int iSub = 0; iSub < MultiDomainDynam::solVecInfo().numDom; ++iSub) {
    reducedInfo.domLen[iSub] = (iSub==0) ? normalizedBasis_.numVec() : 0;
    totLen += reducedInfo.domLen[iSub];
  }

  reducedInfo.len = totLen;
  reducedInfo.setMasterFlag();
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

  normalizedBasis_.projectDown( *d_n, _d_n);
  normalizedBasis_.projectDown( *v_n, _v_n);
  normalizedBasis_.projectDown( *a_n, _a_n);
  normalizedBasis_.projectDown( *v_p, _v_p);
}

void 
DistrExplicitPodProjectionNonLinDynamicBase::updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n1) {
  //update geomState for Fint, but no need to update displacment vector from geometry 

  DistrVector temp1(solVecInfo());
  temp1 = dt_n_h*v_n_h;

  normalizedBasis_.projectUp( temp1, *d_n); 

  geomState->update(*d_n, 2);

  d_n1 += temp1;  //we save the increment vectors for postprocessing

  if(haveRot) { // currently we only need to project the velocity up when there are rotation dofs
                // int the future, there may be other cases in which this is also necessary, e.g. viscoelastic materials
    normalizedBasis_.projectUp(v_n_h, *v_n);
    geomState->setVelocity(*v_n, 2);
  }

}

void DistrExplicitPodProjectionNonLinDynamicBase::getConstForce(DistrVector& v)
{
  int nr = domain->nNeumannModal();
  if(nr) {
    filePrint(stderr, " ... Using Reduced Constant Force   ...\n");
    BCond* nbcModal = domain->getNBCModal();
    v.zero();
    for(int i=0; i<nr; ++i) v[nbcModal[i].nnum] = nbcModal[i].val;
  }
  else {
    //we really don't need to project down here since cnst_fBig is stored inside the probDesc class
    //just a formality. 
    MultiDomainDynam::getConstForce(*cnst_fBig);
    normalizedBasis_.projectDown(*cnst_fBig,v);
  }
  cnst_fBig->zero();
}

void
DistrExplicitPodProjectionNonLinDynamicBase::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {
  //Build internal force and project into reduced coordinates

  MultiDomainDynam::getInternalForce( *d_n, *fInt, t, tIndex);
  //compute residual here to prevent having to project into reduced basis twice
  *a_n = *fInt - *fExt;

  if(haveRot) {
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

  std::auto_ptr<DistrGalerkinProjectionSolver> solver(new DistrGalerkinProjectionSolver(normalizedBasis_));

  haveRot = geomState->getHaveRot();
  if(haveRot) {
    fullMassSolver = result->dynMat;
  }
  else {
    delete result->dynMat;
  }

  result->dynMat = solver.release();

  return result;
}

} // end namespace Rom
