#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>

#include <Threads.d/PHelper.h>
#include <Feti.d/Feti.h>
#include <Timers.d/GetTime.h>
#include <Paral.d/DomainGroupTask.h>
#include <Utils.d/Memory.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/pstress.h>
#include <Driver.d/DynamProbType.h>
#include <Driver.d/SubDomain.h>
#include <Paral.d/MDDynam.h>
#include <Driver.d/GeoSource.h>
#include <Paral.d/SubDOp.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>

extern const char* problemTypeMessage[];

extern FILE *debugFile;

//#define SERIALIZED_OUTPUT

template<class Scalar>
GenDistrDomain<Scalar>::GenDistrDomain(Domain *d) : GenDecDomain<Scalar>(d)
{
  initialize();
  this->myCPU = this->communicator->cpuNum();
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::initialize()
{
  numRes = 0; 
  x = 0;
  masterFlag = 0; 
  numFlags = 0; 
  nodeOffsets = 0; 
  elemNodeOffsets = 0; 
  nodePat = 0;
  masterStress = 0;
  elemOffsets = 0;
}

template<class Scalar>
GenDistrDomain<Scalar>::~GenDistrDomain()
{
  if(masterFlag) { 
    for(int i=0; i<this->numSub; ++i)
      if(masterFlag[i]) { delete [] masterFlag[i]; masterFlag[i] = 0; }
    delete [] masterFlag; masterFlag = 0; 
  }
  if(numFlags) { delete [] numFlags; numFlags = 0; }
  if(nodeOffsets) { delete [] nodeOffsets; nodeOffsets = 0; }
  if(elemNodeOffsets) { delete [] elemNodeOffsets; elemNodeOffsets = 0; }
  if(numRes) { delete [] numRes; numRes = 0; }
  if(masterStress) { delete masterStress; masterStress = 0;} 
  if(nodePat) { delete nodePat; nodePat = 0; }
  if(elemOffsets) { delete [] elemOffsets; elemOffsets = 0; }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::initPostPro()
{
  //geoSource->setNumNodalOutput();

  if(geoSource->getNumOutInfo()) {
#ifdef DISTRIBUTED
/* moved to DecDomain::preProcess
    if(geoSource->getNumNodalOutput()) {
      for(int i=0; i<this->numSub; ++i)
        geoSource->distributeOutputNodesX(this->subDomain[i], this->nodeToSub);
    }
*/
    createMasterFlag();
    createOutputOffsets();
    makeMasterInfo();
#endif

/*`
    if(!this->elemToNode) { // check if elemToNode is required
      OutputInfo *oinfo = geoSource->getOutputInfo();
      for(int iInfo = 0; iInfo < geoSource->getNumOutInfo(); iInfo++) {
        if(oinfo[iInfo].averageFlg == 0) {
          this->createElemToNode();
          break;
        }
      }
    }
*/ 
    if(!this->elemToNode && geoSource->elemOutput()) this->createElemToNode();
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::makeMasterInfo()
{
  masterInfo.domLen = new int[this->numSub];
  masterInfo.numDom = this->numSub;
  for(int iSub = 0; iSub < this->numSub; ++iSub) {
    masterInfo.domLen[iSub] = numFlags[iSub];
  }
#ifdef DISTRIBUTED
  masterInfo.computeOffsets();
#else
  masterInfo.setMasterFlag();
#endif
}

template<class Scalar>
void
GenDistrDomain<Scalar>::postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f,
                                       double eigV, GenDistrVector<Scalar> *aeroF, int x,
                                       GenMDDynamMat<Scalar> *dynOps, SysState<GenDistrVector<Scalar> > *distState, int ndflag)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;
  int iOut_ffp = -1;

  int outLimit = geoSource->getOutLimit();
  if(numOutInfo && x == 0 && ndflag == 0 && !domain->solInfo().isDynam())
    filePrint(stderr," ... Postprocessing                 ...\n");
  if(!masterFlag) initPostPro();

  int iSub;

  // initialize and merge displacements from subdomains into cpu array
  DistSVec<Scalar, 11> disps(this->nodeInfo);
  DistSVec<Scalar, 11> masterDisps(masterInfo);
  disps = 0;
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
    Scalar *bcx = this->subDomain[iSub]->getBcx();
    this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(xyz, u.subData(iSub), bcx);
  }
  if(domain->solInfo().isCoupled && domain->solInfo().isMatching) unify(disps); // PJSA 1-17-08 make sure master has both fluid and structure solutions before reducing
  disps.reduce(masterDisps, masterFlag, numFlags);

  // initialize and merge aeroelastic forces
  DistSVec<Scalar, 6> aerof(this->nodeInfo);
  DistSVec<Scalar, 6> masterAeroF(masterInfo);
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    GenDistrVector<Scalar> assembledAeroF(*aeroF);
    this->ba->assemble(assembledAeroF);
    aerof = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedAeroF, assembledAeroF.subData(iSub));
    }
    aerof.reduce(masterAeroF, masterFlag, numFlags);
  }

  // initialize and merge velocities & accelerations
  DistSVec<Scalar, 11> vels(this->nodeInfo), accs(this->nodeInfo);
  DistSVec<Scalar, 11> masterVels(masterInfo), masterAccs(masterInfo);
  if(distState) {
    GenDistrVector<Scalar> *v_n = &distState->getVeloc();
    GenDistrVector<Scalar> *a_n = &distState->getAccel();
    vels = 0; accs = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedVel)[11] = (Scalar (*)[11]) vels.subData(iSub);
      Scalar (*mergedAcc)[11] = (Scalar (*)[11]) accs.subData(iSub);
      double *vcx = this->subDomain[iSub]->getVcx();
      Scalar *vcx_scalar = (Scalar *) dbg_alloca(this->subDomain[iSub]->numdof()*sizeof(Scalar));
      for(int i=0; i<this->subDomain[iSub]->numdof(); ++i) vcx_scalar[i] = vcx[i];
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedVel, v_n->subData(iSub), vcx_scalar);
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedAcc, a_n->subData(iSub));
    }
    vels.reduce(masterVels, masterFlag, numFlags);
    accs.reduce(masterAccs, masterFlag, numFlags);
  }

  // compute current time (or frequency in the case of a helmholtz problem)
  double time;
  if(geoSource->isShifted() && domain->probType() != SolverInfo::Modal 
                            && domain->probType() != SolverInfo::Dynamic) {
    time = domain->getFrequencyOrWavenumber();
    if(domain->solInfo().doFreqSweep) x = this->outFreqCount++;
  } else if(domain->probType() == SolverInfo::Modal) {
    time = eigV;
    if(domain->solInfo().doEigSweep) x = this->outEigCount++;
  }
  else time = x*domain->solInfo().getTimeStep();

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

// RT - serialize the OUTPUT,  PJSA - stress output doesn't work with serialized output. need to reconsider
#ifdef SERIALIZED_OUTPUT
for(int iCPU = 0; iCPU < this->communicator->size(); iCPU++) {
 this->communicator->sync();
 if(this->communicator->cpuNum() == iCPU) {
#endif

  // open binary output files
  if(x==0) {
    if(!numRes) numRes = new int[numOutInfo];
    for(int i=0; i<numOutInfo; ++i) numRes[i] = 0;
  }

  if((x==0) || (outLimit > 0 && x%outLimit == 0)) { // PJSA 3-31-06
#ifdef DISTRIBUTED

    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].type == OutputInfo::Farfield || oinfo[iInfo].type == OutputInfo::AeroForce) { 
        int oI = iInfo;
        if(this->firstOutput) { geoSource->openOutputFiles(0,&oI,1); } 
        continue;
      }
      else if(oinfo[iInfo].nodeNumber == -1 && this->firstOutput) { // PJSA only need to call this the first time
        if(this->communicator->cpuNum() == 0) geoSource->createBinaryOutputFile(iInfo,this->localSubToGl[0],x);
        else geoSource->setHeaderLen(iInfo);
      }
    }
#ifndef SERIALIZED_OUTPUT
    this->communicator->sync();
#endif

    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].nodeNumber == -1 && oinfo[iInfo].type != OutputInfo::Farfield && oinfo[iInfo].type != OutputInfo::AeroForce) {
        numRes[iInfo] = 0;
        for(iSub = 0; iSub < this->numSub; iSub++) {
          int glSub = this->localSubToGl[iSub];
          if(oinfo[iInfo].dataType == 1) {
            geoSource->outputRange(iInfo, masterFlag[iSub],
                       this->subDomain[iSub]->numNode(), glSub, nodeOffsets[iSub], x);
          }
          else {
            geoSource->outputRange(iInfo, this->subDomain[iSub]->getGlElems(),
                       this->subDomain[iSub]->numElements(), glSub,
                       elemOffsets[iSub], x);
          }
        }
      }
    }


    if(x==0) { // always put single node output in one file regardless of outLimit
      for(iSub = 0; iSub < this->numSub; iSub++)  {
        int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
        if(nOutNodes) {
          if(this->firstOutput) geoSource->openOutputFiles(this->subDomain[iSub]->getOutputNodes(), 
                   this->subDomain[iSub]->getOutIndex(), nOutNodes);                    
        }
      }
    }

#else
    if(x == 0 && this->firstOutput) geoSource->openOutputFiles();  // opens all output files
#endif
  }

  if (!nodePat) makeNodePat();

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(oinfo[iOut].interval == 0 || x % oinfo[iOut].interval != 0)
      continue;

    // update number of results
    numRes[iOut]++; 

    if(oinfo[iOut].ndtype != ndflag) continue;
    if(ndflag !=0 && oinfo[iOut].type != OutputInfo::Disp6DOF && oinfo[iOut].type !=  OutputInfo::Displacement) continue;

    switch(oinfo[iOut].type)  {

      case OutputInfo::EigenPair:
      case OutputInfo::FreqRespModes:
      case OutputInfo::Displacement:
        getPrimal(disps, masterDisps, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Velocity:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Acceleration:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Disp6DOF:
        getPrimal(disps, masterDisps, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Velocity6:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Accel6:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Temperature:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 6);
        break;
      case OutputInfo::TemperatureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 6);
        break;
      case OutputInfo::AcousticPressure:
      case OutputInfo::EigenPressure:
      case OutputInfo::HelmholtzModes:
      case OutputInfo::Helmholtz:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 7);
        break;
      case OutputInfo::PressureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 7);
        break;
      case OutputInfo::PressureSecondTimeDerivative:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 1, 7);
        break;
      case OutputInfo::StressXX:
        getStressStrain(u, time, x, iOut, SXX);
        break;
      case OutputInfo::StressYY:
        getStressStrain(u, time, x, iOut, SYY);
        break;
      case OutputInfo::StressZZ:
        getStressStrain(u, time, x, iOut, SZZ);
        break;
      case OutputInfo::StressXY:
        getStressStrain(u, time, x, iOut, SXY);
        break;
      case OutputInfo::StressYZ:
        getStressStrain(u, time, x, iOut, SYZ);
        break;
      case OutputInfo::StressXZ:
        getStressStrain(u, time, x, iOut, SXZ);
        break;
      case OutputInfo::StrainXX:
        getStressStrain(u, time, x, iOut, EXX);
        break;
      case OutputInfo::StrainYY:
        getStressStrain(u, time, x, iOut, EYY);
        break;
      case OutputInfo::StrainZZ:
        getStressStrain(u, time, x, iOut, EZZ);
        break;
      case OutputInfo::StrainXY:
        getStressStrain(u, time, x, iOut, EXY);
        break;
      case OutputInfo::StrainYZ:
        getStressStrain(u, time, x, iOut, EYZ);
        break;
      case OutputInfo::StrainXZ:
        getStressStrain(u, time, x, iOut, EXZ);
        break;
      case OutputInfo::StressVM:
        getStressStrain(u, time, x, iOut, VON);
        break;
      case OutputInfo::StrainVM:
        getStressStrain(u, time, x, iOut, STRAINVON);
        break;
      case OutputInfo::ContactPressure:
        getStressStrain(u, time, x, iOut, CONPRESS);
        break;
      case OutputInfo::Damage:
        getStressStrain(u, time, x, iOut, DAMAGE);
        break;
      case OutputInfo::EffPStrn:
        getStressStrain(u, time, x, iOut, EFFPSTRN);
        break;
      case OutputInfo::HardVar:
        getStressStrain(u, time, x, iOut, HARDVAR);
        break;
      case OutputInfo::StressPR1:
        getPrincipalStress(u, time, x, iOut, PSTRESS1);
        break;
      case OutputInfo::StressPR2:
        getPrincipalStress(u, time, x, iOut, PSTRESS2);
        break;
      case OutputInfo::StressPR3:
        getPrincipalStress(u, time, x, iOut, PSTRESS3);
        break;
      case OutputInfo::StrainPR1:
        getPrincipalStress(u, time, x, iOut, PSTRAIN1);
        break;
      case OutputInfo::StrainPR2:
        getPrincipalStress(u, time, x, iOut, PSTRAIN2);
        break;
      case OutputInfo::StrainPR3:
        getPrincipalStress(u, time, x, iOut, PSTRAIN3);
        break;
      case OutputInfo::InXForce:
        getElementForce(u, time, x, iOut, INX);
        break;
      case OutputInfo::InYForce:
        getElementForce(u, time, x, iOut, INY);
        break;
      case OutputInfo::InZForce:
        getElementForce(u, time, x, iOut, INZ);
        break;
      case OutputInfo::AXMoment:
        getElementForce(u, time, x, iOut, AXM);
        break;
      case OutputInfo::AYMoment:
        getElementForce(u, time, x, iOut, AYM);
        break;
      case OutputInfo::AZMoment:
        getElementForce(u, time, x, iOut, AZM);
        break;
      case OutputInfo::DispX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 0);
        break;
      case OutputInfo::DispY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 1);
        break;
      case OutputInfo::DispZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 2);
        break;
      case OutputInfo::RotX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 3);
        break;
      case OutputInfo::RotY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 4);
        break;
      case OutputInfo::RotZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 5);
        break;
      case OutputInfo::DispMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);//DofSet::max_known_nonL_dof
          Scalar *dispMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            dispMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                               xyz[iNode][1]*xyz[iNode][1] +
                                               xyz[iNode][2]*xyz[iNode][2]);
          }
          geoSource->writeNodeScalarToFile(dispMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] dispMod;
        }
        break;
      case OutputInfo::RotMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);//DofSet::max_known_nonL_dof
          Scalar *rotMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            rotMod[iNode] = ScalarTypes::sqrt(xyz[iNode][3]*xyz[iNode][3] +
                                              xyz[iNode][4]*xyz[iNode][4] +
                                              xyz[iNode][5]*xyz[iNode][5]);
          }
          geoSource->writeNodeScalarToFile(rotMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] rotMod;
        }
        break;
      case OutputInfo::TotMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);//DofSet::max_known_nonL_dof
          Scalar *totMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            totMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                              xyz[iNode][1]*xyz[iNode][1] +
                                              xyz[iNode][2]*xyz[iNode][2] +
                                              xyz[iNode][3]*xyz[iNode][3] +
                                              xyz[iNode][4]*xyz[iNode][4] +
                                              xyz[iNode][5]*xyz[iNode][5]);
          }
          geoSource->writeNodeScalarToFile(totMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] totMod;
        }
        break;
      case OutputInfo::Farfield: 
        domain->nffp = oinfo[iOut].interval;
        iOut_ffp = iOut; // PJSA 3-1-2007 buildFFP doesn't work with serialized output
        //this->buildFFP(u,oinfo[iOut].filptr);
        break;
      case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
      case OutputInfo::AeroXForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 0);
        break;
      case OutputInfo::AeroYForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 1);
        break;
      case OutputInfo::AeroZForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 2);
        break;
      case OutputInfo::AeroXMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 3);
        break;
      case OutputInfo::AeroYMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 4);
        break;
      case OutputInfo::AeroZMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 5);
        break;
      case OutputInfo::EigenSlosh:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 10);
        break;
      case OutputInfo::YModulus:
        this->getElementAttr(iOut,YOUNG, time);
        break;
      case OutputInfo::MDensity:
        this->getElementAttr(iOut,MDENS, time);
        break;
      case OutputInfo::Thicknes:
        this->getElementAttr(iOut,THICK, time);
        break;
      case OutputInfo::TDEnforcement: {
        DistSVec<double, 1> all_data(this->nodeInfo); all_data = 0;
        double **sub_data = new double * [this->numSub];
        for(iSub = 0; iSub < this->numSub; ++iSub) sub_data[iSub] = (double *) all_data.subData(iSub);
        for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
          domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[iOut].tdenforc_var, sub_data, this->numSub, this->subDomain);
        }
        DistSVec<double, 1> master_data(masterInfo);
        all_data.reduce(master_data, masterFlag, numFlags);
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          geoSource->writeNodeScalarToFile((double *) master_data.subData(iSub), master_data.subSize(iSub), this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
        }
        delete [] sub_data;
      } break;
      default:
        filePrint(stderr," *** WARNING: Output case %d not implemented \n", iOut);
        break;
    }
  }
  this->firstOutput = false;
// RT - serialize the OUTPUT
#ifdef SERIALIZED_OUTPUT
  }
}
this->communicator->sync();
#endif
  if(iOut_ffp > -1) this->buildFFP(u,oinfo[iOut_ffp].filptr); // PJSA 3-1-2007 buildFFP doesn't work with serialized output
}


template<class Scalar>
void
GenDistrDomain<Scalar>::getPrimal(DistSVec<Scalar, 11> &disps, DistSVec<Scalar, 11> &masterDisps, 
                                  double time, int x, int fileNumber, int ndof, int startdof)
{
  // this function outputs the primal variables: displacement, temperature, pressure
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  int iSub;

  if(oinfo.nodeNumber == -1) { // output binary data for all nodes
    for(iSub = 0; iSub < this->numSub; ++iSub)   
      geoSource->writeNodeVectorToFile(masterDisps(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber],
                                       time, ndof, startdof, masterFlag[iSub]); 
  }
  else {  
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes) {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++) {
          if(outIndex[iNode] == fileNumber) {
            Scalar (*nodeDisp)[11] = (Scalar (*)[11]) disps.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            switch(ndof) {
              case 6:
                geoSource->outputNodeVectors6(fileNumber, nodeDisp + outNodes[iNode], 1, time);
                break;
              case 3:
                geoSource->outputNodeVectors(fileNumber, nodeDisp + outNodes[iNode], 1, time);
                break;
              case 1:
                geoSource->outputNodeScalars(fileNumber, nodeDisp[outNodes[iNode]]+startdof, 1, time);
                break;
              default:
                filePrint(stderr, " *** WARNING: single node primal output not supported for ndof = %d \n", ndof);
                break;
            }
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getAeroForceScalar(DistSVec<Scalar, 6> &aerof, DistSVec<Scalar, 6> &masterAeroF,
                                           double time, int x, int fileNumber, int dof)
{
  // this function outputs an aeroelastic force scalar  
  // dof == 0 is X Force, dof == 1 is Y Force,  etc.
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  int iSub;

  if(oinfo.nodeNumber == -1) { // output binary data for all nodes
    for(iSub = 0; iSub < this->numSub; ++iSub)
      geoSource->writeNodeVectorToFile(masterAeroF(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber],
                                       time, 1, dof, masterFlag[iSub]);
  }
  else { // output ascii data at selected node or nodes only
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes)  {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++)
          if(outIndex[iNode] == fileNumber)  {
            Scalar (*nodeAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            geoSource->outputNodeScalars(fileNumber, nodeAeroF[outNodes[iNode]]+dof, 1, time);
          }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1) this->sizeSfemStress = masterInfo.totLen(); 
  else if(avgnum == 0) {  // element-based output
  this->sizeSfemStress = 0;
/*   Connectivity *elemToNode = new Connectivity(domain->getEset());
   int numele = geoSource->getNumAttributes();  // number of elements; another option domain->numElements();
   for(int iele=0; iele<numele; ++iele)   {
//     cerr << "number of nodes in this element  = " << elemToNode->num(iele) << endl;
     sizeSfemStress = sizeSfemStress + elemToNode->num(iele); // add number of nodes for each element
   }*/
  }
  else {
   cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << endl;
   this->sizeSfemStress = 0;
  }
}

template<class Scalar>
Scalar*
GenDistrDomain<Scalar>::getSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1) return masterStress->data(); // node-based
  else if(avgnum == 0) return 0; // element-based YYY DG Implement
  else  { cerr << "avgnum = " << avgnum << " not implemented in Domain::getSfemStress()" << endl; return 0; }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::updateSfemStress(Scalar* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  //int numNodes = masterInfo.totLen(); //  size(masterStress)

  if(avgnum == 1)  masterStress->setNewData(str); // YYY DG
  else if(avgnum == 0) cerr << "updateSfemStress for element not yet implemented" << endl; // YYY DG
  else {cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << endl;}
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getStressStrain(GenDistrVector<Scalar> &u, double time,
                                        int x, int fileNumber, int Findex, int printFlag) 
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementStressStrain(u, time, x, fileNumber, Findex, printFlag);
    return;
  }

  DistVec<Scalar> stress(this->nodeInfo);
  DistVec<Scalar> weight(this->nodeInfo);

  stress = 0;
  weight = 0;


  int iSub;
  
  if(printFlag != 2) {
    // each subdomain computes its stress vector
    for (iSub = 0; iSub < this->numSub; ++iSub) {
      this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub), 
  		Findex, stress.subData(iSub), weight.subData(iSub));
    }

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, nodePat, &stress);
    nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, nodePat, &stress);

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, nodePat, &weight);
    nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, nodePat, &weight);

    // Divide stress by weight
    for (iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = stress(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress.subSize(iSub); ++i)
        if(locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
    }

    // reduce the stress vector to just the master quantities
    if(masterStress == 0) masterStress = new DistVec<Scalar>(masterInfo);
    stress.reduce(*masterStress, masterFlag, numFlags);
  }
  

  if(printFlag != 1) {
    // write to file
    for (iSub = 0; iSub < this->numSub; ++iSub) {
      geoSource->writeNodeScalarToFile(masterStress->subData(iSub), masterStress->subSize(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]);
    }
  }

}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementStressStrain(GenDistrVector<Scalar> &u, double time, 
                                               int iter, int fileNumber, int Findex, int printFlag) 
{
  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemStress = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub), Findex, elemStress);
    geoSource->writeElemScalarToFile(elemStress, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemStress;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementPrincipalStress(GenDistrVector<Scalar> &u, double time,
                                                  int iter, int fileNumber, int strIndex)
{
  // PJSA: 3-23-05
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
                                                                                                                                                             
  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++) {  // this could be parallelized
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar **elemAllStress = new Scalar * [6];
    for(i=0; i<6; ++i) elemAllStress[i] = new Scalar[numElemNodes];
                                                                                                                               
    // Compute Each Required Stress (all 6) using same routines as for
    // individual stresses
    int str_loop;
    int Findex;
    for(str_loop = 0; str_loop < 6; ++str_loop) {
      // get current stress/strain index
      Findex = strDir[str_loop];

      // compute stress vector
      this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub),
                                           Findex, elemAllStress[str_loop]);
    }

    // ... CALCULATE PRINCIPALS AT EACH NODE
    Scalar svec[6], pvec[3];
    Scalar *elemPVec = new Scalar[numElemNodes];
    for(i = 0; i < numElemNodes; ++i) {
      for(j = 0; j < 6; ++j)
        svec[j] = elemAllStress[j][i];
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec, pvec);
      elemPVec[i] = pvec[strInd];
    }

    geoSource->writeElemScalarToFile(elemPVec, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());

    delete [] elemPVec;
    for(i=0; i<6; ++i) delete [] elemAllStress[i];
    delete [] elemAllStress;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                                                  int x, int fileNumber, int strIndex)
{
  // PJSA: 3-23-05 Non-linear version
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++) {  // this could be parallelized
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar **elemAllStress = new Scalar * [6];
    for(i=0; i<6; ++i) elemAllStress[i] = new Scalar[numElemNodes];
                                                                                                                                                           
    // Compute Each Required Stress (all 6) using same routines as for
    // individual stresses
    int str_loop;
    int Findex;
    for(str_loop = 0; str_loop < 6; ++str_loop) {
      // get current stress/strain index
      Findex = strDir[str_loop];
                                                                                                                                                           
      // compute stress vector
      this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber,
                                           Findex, elemAllStress[str_loop]);
    }
    // ... CALCULATE PRINCIPALS AT EACH NODE
    Scalar svec[6], pvec[3];
    Scalar *elemPVec = new Scalar[numElemNodes];
    for(i = 0; i < numElemNodes; ++i) {
      for(j = 0; j < 6; ++j)
        svec[j] = elemAllStress[j][i];
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec, pvec);
      elemPVec[i] = pvec[strInd];
    }

    geoSource->writeElemScalarToFile(elemPVec, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, x,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());

    delete [] elemPVec;
    for(i=0; i<6; ++i) delete [] elemAllStress[i];
    delete [] elemAllStress;
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getPrincipalStress(GenDistrVector<Scalar> &u, double time, int x, 
                                           int fileNumber, int strIndex)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(u, time, x, fileNumber, strIndex);
    return;
  }

  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress
  DistVec<Scalar> **stress = new DistVec<Scalar>*[6];
  DistVec<Scalar> weight(this->nodeInfo);

  int iSub;
  int str_loop;
  for(str_loop = 0; str_loop < 6; ++str_loop)
    stress[str_loop] = new DistVec<Scalar> (this->nodeInfo);  

  // each subdomain computes its stress/strain vector

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int Findex;

  for(str_loop = 0; str_loop < 6; ++str_loop) {

    // get current stress/strain index
    Findex = strDir[str_loop];

    // Initialize distributed vector to zero
    *stress[str_loop] = 0;
    weight = 0;

    for (iSub = 0; iSub < this->numSub; ++iSub) {
      this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub),
                Findex, stress[str_loop]->subData(iSub), weight.subData(iSub));
    }

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>,
               this->nodePat, stress[str_loop]);
    this->nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, stress[str_loop]);

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
    this->nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

    // Divide stress by weight
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = (*stress[str_loop])(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress[str_loop]->subSize(iSub); ++i)  {
        if (locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
      }
    }
  }

  // Calculate Principals at each node
  Scalar svec[6], pvec[3];
  DistVec<Scalar> allPVec(this->nodeInfo);
  for(iSub = 0; iSub < this->numSub; ++iSub)  {

    Vec<Scalar> &locPVec = allPVec(iSub);
    for(i = 0; i < this->subDomain[iSub]->numNode(); ++i) {

      for(j = 0; j < 6; ++j)
        svec[j] = stress[j]->subData(iSub)[i];
 
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      locPVec[i] = pvec[strInd];
    }
  }

  // reduce stress vector to master quantities
  DistVec<Scalar> masterPVec(masterInfo);
  allPVec.reduce(masterPVec, masterFlag, numFlags);

  // print out data
  for (iSub = 0; iSub < this->numSub; ++iSub)  {
    geoSource->writeNodeScalarToFile(masterPVec.subData(iSub), masterPVec.subSize(iSub), this->localSubToGl[iSub],
                                     nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getElementForce(GenDistrVector<Scalar> &u, double time, int x, 
                                        int fileNumber, int Findex )  
{
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemForce = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeElementForce(u.subData(iSub),
                                         Findex, elemForce);
    geoSource->writeElemScalarToFile(elemForce, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, x,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemForce;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::createMasterFlag() 
{
  // allocate size of master flag
  masterFlag = new int *[this->numSub];
  numFlags = new int[this->numSub];

  int iSub;

  for(iSub = 0; iSub < this->numSub; iSub++)  {
    int numNodes = this->subDomain[iSub]->numNode();
    masterFlag[iSub] = new int[numNodes];
    numFlags[iSub] = numNodes;

    // initialize master flag for each subdomain
    int *clusNodes = this->subDomain[iSub]->getGlNodes();
    for (int iNode = 0; iNode < numNodes; iNode++)
      masterFlag[iSub][iNode] = clusNodes[iNode];
  }

  Connectivity *subToClus = geoSource->getSubToClus();
  for(iSub = 0; iSub < this->numSub; iSub++) {
    // get connected subdomain information
    int *connDoms = this->subDomain[iSub]->getSComm()->subNums;
    int numConns = this->subDomain[iSub]->getSComm()->numNeighb;

    int thisGlSub = this->subDomain[iSub]->subNum();
    int thisCluster = (subToClus) ? (*subToClus)[thisGlSub][0] : 0;

    // loop over connected subdomains in this cluster
    for(int jSub = 0; jSub < numConns; jSub++) {
      // set master flag if connected subdomain global number is smaller
      if(connDoms[jSub] < thisGlSub && (subToClus == 0 || (*subToClus)[connDoms[jSub]][0] == thisCluster)) {

        // get shared nodes
        Connectivity *sharedNodes = this->subDomain[iSub]->getSComm()->sharedNodes;
        for(int iNode = 0; iNode < sharedNodes->num(jSub); iNode++) {
          int shNode = (*sharedNodes)[jSub][iNode];
          if(masterFlag[iSub][shNode] >= 0) {
            masterFlag[iSub][shNode] = -1;
            numFlags[iSub]--;
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::createOutputOffsets() 
{
  Connectivity *clusToSub = geoSource->getClusToSub();
  Connectivity *subToClus = geoSource->getSubToClus();
  nodeOffsets = new int[this->numSub];
  elemNodeOffsets = new int[this->numSub];
  elemOffsets = new int[this->numSub];

  int iSub, iCluster;
  int numClusters = geoSource->getNumClusters();
  if(numClusters == 0) return;
  // create cluster array of offsets
  int **clNodeOffsets = new int *[numClusters];
  int **clElemNodeOffsets = new int *[numClusters];
  int **clElemOffsets = new int *[numClusters];  // PJSA

  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    clNodeOffsets[iCluster] = new int[clusToSub->num(iCluster)];
    clElemNodeOffsets[iCluster] = new int[clusToSub->num(iCluster)];
    clElemOffsets[iCluster] = new int[clusToSub->num(iCluster)];
    for(int j=0; j<clusToSub->num(iCluster); ++j) { // PJSA: initialize to zero so global sum will work
      clNodeOffsets[iCluster][j] = 0;
      clElemNodeOffsets[iCluster][j] = 0;
      clElemOffsets[iCluster][j] = 0;
    }
  }

  // populate glOffsets w/this mpi's data
  for(iSub = 0; iSub < this->numSub; iSub++) {
    int glSub = this->subDomain[iSub]->subNum();
    int clusNum = (*subToClus)[glSub][0];
    for(int jSub = 0; jSub < clusToSub->num(clusNum); jSub++)
      if(glSub == (*clusToSub)[clusNum][jSub]) {
        clNodeOffsets[clusNum][jSub] = numFlags[iSub];
        clElemNodeOffsets[clusNum][jSub] = this->subDomain[iSub]->countElemNodes();
        clElemOffsets[clusNum][jSub] = this->subDomain[iSub]->numElements();
      }
  }

  // sum up all offsets in all mpi processes
  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    this->communicator->globalSum(clusToSub->num(iCluster), clNodeOffsets[iCluster]);
    this->communicator->globalSum(clusToSub->num(iCluster), clElemNodeOffsets[iCluster]);
    this->communicator->globalSum(clusToSub->num(iCluster), clElemOffsets[iCluster]);
  }

  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    int nOffset = 0;
    int enOffset = 0;
    int eOffset = 0;
    for(iSub = 0; iSub < clusToSub->num(iCluster); iSub++) {
      int tmpNOff = clNodeOffsets[iCluster][iSub];
      int tmpENOff = clElemNodeOffsets[iCluster][iSub];
      int tmpEOff = clElemOffsets[iCluster][iSub];
      clNodeOffsets[iCluster][iSub] = nOffset;
      clElemNodeOffsets[iCluster][iSub] = enOffset;
      clElemOffsets[iCluster][iSub] = eOffset;
      nOffset += tmpNOff;
      enOffset += tmpENOff;
      eOffset += tmpEOff;
    }
  }
   
  for(iSub = 0; iSub < this->numSub; iSub++) {
    int glSub = this->subDomain[iSub]->subNum();
    int clusNum = (*subToClus)[glSub][0];
    for(int jSub = 0; jSub < clusToSub->num(clusNum); jSub++)
      if(glSub == (*clusToSub)[clusNum][jSub]) {
        nodeOffsets[iSub] = clNodeOffsets[clusNum][jSub];
        elemNodeOffsets[iSub] = clElemNodeOffsets[clusNum][jSub];
        elemOffsets[iSub] = clElemOffsets[clusNum][jSub];
      }
  }

  if(clNodeOffsets) { for(int i=0; i<numClusters; i++) delete [] clNodeOffsets[i]; delete [] clNodeOffsets; }
  if(clElemNodeOffsets) { for(int i=0; i<numClusters; i++) delete [] clElemNodeOffsets[i]; delete [] clElemNodeOffsets; }
  if(clElemOffsets) { for(int i=0; i<numClusters; i++) delete [] clElemOffsets[i]; delete [] clElemOffsets; }
}

// ---------------------------------------------------------------
// Nonlinear DistrDomain functions
// ----------------------------------------------------------------

template<class Scalar>
void
GenDistrDomain<Scalar>::postProcessing(DistrGeomState *geomState, Corotator ***allCorot, double time, SysState<GenDistrVector<Scalar> > *distState,
                                       GenDistrVector<Scalar> *aeroF)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  int outLimit = geoSource->getOutLimit();
  if(numOutInfo && x == 0 && !domain->solInfo().isDynam())
    filePrint(stderr," ... Postprocessing                 ...\n");
  if(!masterFlag) initPostPro();

  int iSub;

  // initialize and merge displacements from subdomains into cpu array
  DistSVec<Scalar, 11> disps(this->nodeInfo);
  DistSVec<Scalar, 11> masterDisps(masterInfo);
  disps = 0;
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);//DofSet::max_known_nonL_dof
    this->subDomain[iSub]->mergeDistributedNLDisp(xyz, (*geomState)[iSub]);
  }
  if(domain->solInfo().isCoupled && domain->solInfo().isMatching) unify(disps); // PJSA 1-17-08 make sure master has both fluid and structure solutions before reducing
  disps.reduce(masterDisps, masterFlag, numFlags);

  // initialize and merge aeroelastic forces
  DistSVec<Scalar, 6> aerof(this->nodeInfo);
  DistSVec<Scalar, 6> masterAeroF(masterInfo);
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    GenDistrVector<Scalar> assembledAeroF(*aeroF);
    this->ba->assemble(assembledAeroF);
    aerof = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedAeroF, assembledAeroF.subData(iSub));
    }
    aerof.reduce(masterAeroF, masterFlag, numFlags);
  }

  // initialize and merge velocities & accelerations
  DistSVec<Scalar, 11> vels(this->nodeInfo), accs(this->nodeInfo);
  DistSVec<Scalar, 11> masterVels(masterInfo), masterAccs(masterInfo);
  if(distState) {
    GenDistrVector<Scalar> *v_n = &distState->getVeloc();
    GenDistrVector<Scalar> *a_n = &distState->getAccel();
    vels = 0; accs = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedVel)[11] = (Scalar (*)[11]) vels.subData(iSub);
      Scalar (*mergedAcc)[11] = (Scalar (*)[11]) accs.subData(iSub);
      double *vcx = this->subDomain[iSub]->getVcx();
      Scalar *vcx_scalar = (Scalar *) dbg_alloca(this->subDomain[iSub]->numdof()*sizeof(Scalar));
      for(int i=0; i<this->subDomain[iSub]->numdof(); ++i) vcx_scalar[i] = vcx[i];
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedVel, v_n->subData(iSub), vcx_scalar);
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedAcc, a_n->subData(iSub));
    }
    vels.reduce(masterVels, masterFlag, numFlags);
    accs.reduce(masterAccs, masterFlag, numFlags);
  }

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

// RT - serialize the OUTPUT, PJSA - stress output doesn't work with serialized output. need to reconsider
#ifdef SERIALIZED_OUTPUT
for(int iCPU = 0; iCPU < this->communicator->size(); iCPU++) {
 this->communicator->sync();
 if(this->communicator->cpuNum() == iCPU) {
#endif

  // open binary output files
  if(x==0) {
    if(!numRes) numRes = new int[numOutInfo];
    for(int i=0; i<numOutInfo; ++i) numRes[i] = 0;
  }

  if((x==0) || (outLimit > 0 && x%outLimit == 0)) { // PJSA 3-31-06
#ifdef DISTRIBUTED

    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].type == OutputInfo::Farfield || oinfo[iInfo].type == OutputInfo::AeroForce) {
        int oI = iInfo;
        if(this->firstOutput) { geoSource->openOutputFiles(0,&oI,1); }
        continue;
      }
      else if(oinfo[iInfo].nodeNumber == -1 && this->firstOutput) { // PJSA only need to call this the first time
        if(this->communicator->cpuNum() == 0) geoSource->createBinaryOutputFile(iInfo,this->localSubToGl[0],x);
        else geoSource->setHeaderLen(iInfo);
      }
    }
#ifndef SERIALIZED_OUTPUT
    this->communicator->sync();
#endif

    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].nodeNumber == -1 && oinfo[iInfo].type != OutputInfo::Farfield && oinfo[iInfo].type != OutputInfo::AeroForce) {
        numRes[iInfo] = 0;
        for(iSub = 0; iSub < this->numSub; iSub++) {
          int glSub = this->localSubToGl[iSub];
          if(oinfo[iInfo].dataType == 1) {
            geoSource->outputRange(iInfo, masterFlag[iSub],
                       this->subDomain[iSub]->numNode(), glSub, nodeOffsets[iSub], x);
          }
          else {
            geoSource->outputRange(iInfo, this->subDomain[iSub]->getGlElems(),
                       this->subDomain[iSub]->numElements(), glSub,
                       elemOffsets[iSub], x);
          }
        }
      }
    }


    if(x==0) { // always put single node output in one file regardless of outLimit
      for(iSub = 0; iSub < this->numSub; iSub++)  {
        int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
        if(nOutNodes) {
          if(this->firstOutput) geoSource->openOutputFiles(this->subDomain[iSub]->getOutputNodes(),
                   this->subDomain[iSub]->getOutIndex(), nOutNodes);
        }
      }
    }

#else
    if(x == 0 && this->firstOutput) geoSource->openOutputFiles();  // opens all output files
#endif
  }

  if (!nodePat) makeNodePat();

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(oinfo[iOut].interval == 0 || x % oinfo[iOut].interval != 0)
      continue;

    // update number of results
    numRes[iOut]++;

    switch(oinfo[iOut].type)  {
      
      case OutputInfo::FreqRespModes:
      case OutputInfo::Displacement:
        getPrimal(disps, masterDisps, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Velocity:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Acceleration:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Disp6DOF:
        getPrimal(disps, masterDisps, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Velocity6:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Accel6:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Temperature:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 6);
        break;
      case OutputInfo::TemperatureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 6);
        break;
      case OutputInfo::StressXX:
        getStressStrain(geomState, allCorot, time, x, iOut, SXX);
        break;
      case OutputInfo::StressYY:
        getStressStrain(geomState, allCorot, time, x, iOut, SYY);
        break;
      case OutputInfo::StressZZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SZZ);
        break;
      case OutputInfo::StressXY:
        getStressStrain(geomState, allCorot, time, x, iOut, SXY);
        break;
      case OutputInfo::StressYZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SYZ);
        break;
      case OutputInfo::StressXZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SXZ);
        break;
      case OutputInfo::StrainXX:
        getStressStrain(geomState, allCorot, time, x, iOut, EXX);
        break;
      case OutputInfo::StrainYY:
        getStressStrain(geomState, allCorot, time, x, iOut, EYY);
        break;
      case OutputInfo::StrainZZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EZZ);
        break;
      case OutputInfo::StrainXY:
        getStressStrain(geomState, allCorot, time, x, iOut, EXY);
        break;
      case OutputInfo::StrainYZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EYZ);
        break;
      case OutputInfo::StrainXZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EXZ);
        break;
      case OutputInfo::StressVM:
        getStressStrain(geomState, allCorot, time, x, iOut, VON);
        break;
      case OutputInfo::StrainVM:
        getStressStrain(geomState, allCorot, time, x, iOut, STRAINVON);
        break;
      case OutputInfo::ContactPressure:
        getStressStrain(geomState, allCorot, time, x, iOut, CONPRESS);
        break;
      case OutputInfo::StressPR1:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS1);
        break;
      case OutputInfo::StressPR2:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS2);
        break;
      case OutputInfo::StressPR3:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS3);
        break;
      case OutputInfo::StrainPR1:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN1);
        break;
      case OutputInfo::StrainPR2:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN2);
        break;
      case OutputInfo::StrainPR3:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN3);
        break;
      case OutputInfo::DispX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 0);
        break;
      case OutputInfo::DispY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 1);
        break;
      case OutputInfo::DispZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 2);
        break;
      case OutputInfo::RotX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 3);
        break;
      case OutputInfo::RotY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 4);
        break;
      case OutputInfo::RotZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 5);
        break;
      case OutputInfo::DispMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[8] = (Scalar (*)[8]) masterDisps.subData(iSub);
          Scalar *dispMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            dispMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                               xyz[iNode][1]*xyz[iNode][1] +
                                               xyz[iNode][2]*xyz[iNode][2]);
          }
          geoSource->writeNodeScalarToFile(dispMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] dispMod;
        }
        break;
      case OutputInfo::RotMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[8] = (Scalar (*)[8]) masterDisps.subData(iSub);
          Scalar *rotMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            rotMod[iNode] = ScalarTypes::sqrt(xyz[iNode][3]*xyz[iNode][3] +
                                              xyz[iNode][4]*xyz[iNode][4] +
                                              xyz[iNode][5]*xyz[iNode][5]);
          }
          geoSource->writeNodeScalarToFile(rotMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] rotMod;
        }
        break;
      case OutputInfo::TotMod:
        for(iSub = 0; iSub < this->numSub; ++iSub) {
          int size = masterDisps.subSize(iSub);
          Scalar (*xyz)[8] = (Scalar (*)[8]) masterDisps.subData(iSub);
          Scalar *totMod = new Scalar[size];
          for(int iNode=0; iNode<size; ++iNode) {
            totMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                              xyz[iNode][1]*xyz[iNode][1] +
                                              xyz[iNode][2]*xyz[iNode][2] +
                                              xyz[iNode][3]*xyz[iNode][3] +
                                              xyz[iNode][4]*xyz[iNode][4] +
                                              xyz[iNode][5]*xyz[iNode][5]);
          }
          geoSource->writeNodeScalarToFile(totMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                           iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          delete [] totMod;
        }
        break;
      case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
      case OutputInfo::AeroXForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 0);
        break;
      case OutputInfo::AeroYForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 1);
        break;
      case OutputInfo::AeroZForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 2);
        break;
      case OutputInfo::AeroXMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 3);
        break;
      case OutputInfo::AeroYMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 4);
        break;
      case OutputInfo::AeroZMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 5);
        break;
      default:
        filePrint(stderr," *** WARNING: Output case %d not implemented for non-linear FETI\n", iOut);
        break;
    }
  }
  x++;
#ifdef SERIALIZED_OUTPUT
 }
}
this->communicator->sync();
#endif
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                        int x, int fileNumber, int Findex)
{
  // non-linear version of getStressStrain
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  if(oinfo.averageFlg == 0) {
    getElementStressStrain(gs, allCorot, time, x, fileNumber, Findex);
    return;
  }

  DistVec<Scalar> stress(this->nodeInfo);
  DistVec<Scalar> weight(this->nodeInfo);

  stress = 0;
  weight = 0;

  int iSub;

  // each subdomain computes its stress vector
  for (iSub = 0; iSub < this->numSub; ++iSub) {
    this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber,
                                         Findex, stress.subData(iSub), weight.subData(iSub));
  }

  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &stress);
  this->nodePat->exchange();
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &stress);

  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
  this->nodePat->exchange();
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

  // Divide stress by weight
  for (iSub = 0; iSub < this->numSub; ++iSub)  {
    Vec<Scalar> &locStress = stress(iSub);
    Vec<Scalar> &locWeight = weight(iSub);
    for(int i = 0; i < stress.subSize(iSub); ++i)
      locStress[i] /= locWeight[i];
  }

  // reduce the stress vector to just the master quantities
  DistVec<Scalar> masterStress(masterInfo);
  stress.reduce(masterStress, masterFlag, numFlags);

  // write to file
  for (iSub = 0; iSub < this->numSub; ++iSub) {
    geoSource->writeNodeScalarToFile(masterStress.subData(iSub), masterStress.subSize(iSub), this->localSubToGl[iSub],
                                     nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
  }
}
  
template<class Scalar>
void
GenDistrDomain<Scalar>::getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                               int iter, int fileNumber, int Findex)
{
  // non-linear version of getElementStressStrain
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemStress = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber, Findex, elemStress);
    geoSource->writeElemScalarToFile(elemStress, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemStress;
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,  
                                           int x, int fileNumber, int strIndex)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(gs, allCorot, time, x, fileNumber, strIndex);
    return;
  }

  // non-linear version of getPrincipalStress
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  DistVec<Scalar> **stress = new DistVec<Scalar>*[6];
  DistVec<Scalar> weight(this->nodeInfo);

  int iSub;
  int str_loop;
  for(str_loop = 0; str_loop < 6; ++str_loop)
    stress[str_loop] = new DistVec<Scalar> (this->nodeInfo);  

  // each subdomain computes its stress/strain vector

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int Findex;

  for(str_loop = 0; str_loop < 6; ++str_loop) {

    // get current stress/strain index
    Findex = strDir[str_loop];

    // Initialize distributed vector to zero
    *stress[str_loop] = 0;
    weight = 0;

    for(iSub = 0; iSub < this->numSub; ++iSub) {
      this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber, 
                Findex, stress[str_loop]->subData(iSub), weight.subData(iSub));
    }

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>,
               this->nodePat, stress[str_loop]);
    this->nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, stress[str_loop]);

    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
    this->nodePat->exchange();
    paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

    // Divide stress by weight
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = (*stress[str_loop])(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress[str_loop]->subSize(iSub); ++i)  {
        if(locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
      }
    }
  }

  // Calculate Principals at each node
  Scalar svec[6], pvec[3];
  DistVec<Scalar> allPVec(this->nodeInfo);
  for(iSub = 0; iSub < this->numSub; ++iSub) {

    Vec<Scalar> &locPVec = allPVec(iSub);
    for(i = 0; i < this->subDomain[iSub]->numNode(); ++i) {

      for(j = 0; j < 6; ++j)
        svec[j] = stress[j]->subData(iSub)[i];
 
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      locPVec[i] = pvec[strInd];
    }
  }

  // reduce stress vector to master quantities
  DistVec<Scalar> masterPVec(masterInfo);
  allPVec.reduce(masterPVec, masterFlag, numFlags);

  // print out data
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    geoSource->writeNodeScalarToFile(masterPVec.subData(iSub), masterPVec.subSize(iSub), this->localSubToGl[iSub],
                                     nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::makeNodePat()
{
  nodePat = new FSCommPattern<Scalar>(this->communicator, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSub; ++i) this->subDomain[i]->setNodeCommSize(nodePat);
  nodePat->finalize();
}

template<class Scalar>
void
GenDistrDomain<Scalar>::unify(DistSVec<Scalar, 11> &vec)
{ 
  // XXXX PJSA 1-17-08: make sure that every subdomain sharing a node has the same solution for all the dofs. Actually the master sub is really the only one that needs it.
  FSCommPattern<Scalar> *pat = new FSCommPattern<Scalar>(this->communicator, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSub; ++i) this->subDomain[i]->setNodeCommSize(pat, 11);
  pat->finalize();
  for(int i=0; i<this->numSub; ++i) this->subDomain[i]->sendNode((Scalar (*)[11]) vec.subData(i), pat);
  pat->exchange();
  for(int i=0; i<this->numSub; ++i) this->subDomain[i]->collectNode((Scalar (*)[11]) vec.subData(i), pat);
  delete pat;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenDistrDomain<Scalar>::getElementAttr(int fileNumber,int iAttr, double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 0) 
    {
      assert(0); // not implemented
      return;
    }

  DistVec<double> props(this->nodeInfo);
  DistVec<double> weight(this->nodeInfo);

  // Initialize distributed vector to zero
  props  = 0;
  weight = 0;
  for(int iSub = 0; iSub < this->numSub; ++iSub) 
    { this->getSubDomain(iSub)->mergeElemProps(props.subData(iSub), weight.subData(iSub), iAttr); }
  
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<double>,
	     this->nodePat, &props);
  this->nodePat->exchange();
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<double>, this->nodePat, &props);
  
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<double>, this->nodePat, &weight);
  this->nodePat->exchange();
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::template addNodalData<double>, this->nodePat, &weight);
  
  // average
  for(int iSub = 0; iSub < this->numSub; ++iSub)  
    {
      Vec<double> &locProps =  props (iSub);
      Vec<double> &locWeight = weight(iSub);
      for(int i = 0; i < props.subSize(iSub); ++i)  
	{
	  if (locWeight[i] != 0.0)
	    { locProps[i] /= locWeight[i]; }
	  else
	    { locProps[i] = 0.0; }
	}
    }
  // reduce props vector to master quantities
  DistVec<double> masterVec(masterInfo);
  props.reduce(masterVec, masterFlag, numFlags);
  
  // print out data
  for(int iSub = 0; iSub < this->numSub; ++iSub)  
    {
      geoSource->writeNodeScalarToFile(masterVec.subData(iSub),
				       masterVec.subSize(iSub), this->localSubToGl[iSub],
				       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]);
    }
  return;
}
