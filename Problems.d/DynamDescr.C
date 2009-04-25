#include <stdlib.h>
#include <stdio.h>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>

#include <Problems.d/DynamDescr.h>

#include <Math.d/FullMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Skyline.d/SkyMatrix.h>

#include <Utils.d/dofset.h>
#include <Solvers.d/Solver.h>
#include <Element.d/State.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>

#include <Control.d/ControlInterface.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/BinFileHandler.h>

#include <Utils.d/linkfc.h>
#include <Driver.d/GeoSource.h>

typedef FSFullMatrix FullMatrix;

// ALL COMMENTED FUNCTIONS ARE IN THE FILE DYNAMDESCRTEM.C - JF

//-------------------------------------------------------------------------------------------------------------
template<>
PitaDynamMat *
SingleDomainDynamic<double>::buildPitaOps(double cM, double cC, double cK, double cM_Dt, double cC_Dt, double cK_Dt, bool Coarse)
{
 /* CD: for Pita, factorize amplification matrix with dt and Dt.
        dMat inherits a DynamMat dynMat that contains: K, Kuc, M, 
        Muc, Mcc, C, Cuc and a solver dynMat computed with the 
        time step dt. 
        dMat contains too a solver computed with the time step Dt 
        ie dynMat_Dt.
        All its elements are computed in domain->buildPitaOps
        cf Driver.d/Domain.h and Driver.d/OpMake.C 
  */

 PitaDynamMat *dMat; 

 if (Coarse) {
 
 dMat = new  PitaDynamMat;

 dMat->K = domain->constructDBSparseMatrix<double>();
 dMat->M = domain->constructDBSparseMatrix<double>();
 dMat->Muc = domain->constructCuCSparse<double>();
 dMat->Mcc = domain->constructCCSparse<double>();
 SparseMatrix *Kuc = domain->constructCuCSparse<double>();

 // Rayleigh Damping coefficients
 double alpha = domain->solInfo().alphaDamp;
 double beta  = domain->solInfo().betaDamp;

 // Damping Matrix: C = alpha*M + beta*K
 if(alpha != 0.0 || beta != 0.0) {
    dMat->C   = domain->constructDBSparseMatrix<double>();
    dMat->Cuc = domain->constructCuCSparse<double>();
 }
 
 Rbm *rigidBodyModes = 0;
 
 int useProjector = domain->solInfo().filterFlags;
 int useGrbm = domain->solInfo().rbmflg;

 if (useGrbm || useProjector)
   rigidBodyModes = domain->constructRbm();

 if (!useGrbm || (getTimeIntegration() != 1) )
   domain->buildPitaOps<double>(*dMat, Kuc, cK, cM, cC, cK_Dt, cM_Dt, cC_Dt, 0, kelArray);
 else
   domain->buildPitaOps<double>(*dMat, Kuc, cK, cM, cC, cK_Dt, cM_Dt, cC_Dt, rigidBodyModes, kelArray);

 // Filter RBM

 if(useProjector == 1)
    fprintf(stderr," ... RBMfilter Level 1 Requested    ...\n");
 if(useProjector == 2)
    fprintf(stderr," ... RBMfilter Level 2 Requested    ...\n");

 if(useProjector)
   projector_prep(rigidBodyModes, dMat->M);

 // Modal decomposition preprocessing

 int decompFlag = domain->solInfo().modeDecompFlag;
 if (decompFlag) {
   fprintf(stderr," ... Modal decomposition requested ...\n");
   modeDecompPreProcess(dMat->M);
 }
 
 kuc    = Kuc;
// solver = dMat->dynMat;
 
 }else{
 dMat = (PitaDynamMat *) buildOps(cM, cC, cK); // Unsafe cast !!!
 }

 return dMat;

 fflush(stderr);

}

// SDDynamPostProcessor implementation

SDDynamPostProcessor::SDDynamPostProcessor(Domain *d, double *_bcx, double *_vcx,
                                           StaticTimers *_times)
{ domain = d; bcx = _bcx; vcx = _vcx; times = _times; }

SDDynamPostProcessor::~SDDynamPostProcessor() {
  geoSource->closeOutputFiles();
}

void
SDDynamPostProcessor::openOutputFiles() {
  geoSource->openOutputFiles();
}

void
SDDynamPostProcessor::openOutputFilesForPita(int timeSliceRank) {
  geoSource->openOutputFilesForPita(timeSliceRank);
}

void
SDDynamPostProcessor::closeOutputFiles() {
  geoSource->closeOutputFiles();
}

void
SDDynamPostProcessor::closeOutputFilesForPita(int timeSliceRank) {
  geoSource->closeOutputFilesForPita(timeSliceRank);
}

double 
SDDynamPostProcessor::getKineticEnergy(Vector & vel, SparseMatrix * gMass) {
  return domain->getKineticEnergy(vel,gMass); 
}

void
SDDynamPostProcessor::dynamOutput(int tIndex, DynamMat& dMat, Vector& ext_f, Vector *aeroForce, SysState<Vector> &state) {
  startTimerMemory(times->output, times->memoryOutput);
  
  this->fillBcxVcx(tIndex);

  // PJSA 4-9-08 ext_f passed here may not be for the correct time
  domain->dynamOutput(tIndex, bcx, dMat, ext_f, *aeroForce, state.getDisp(), state.getVeloc(),
                      state.getAccel(), state.getPrevVeloc(), vcx);

  stopTimerMemory(times->output, times->memoryOutput);
}

void
SDDynamPostProcessor::pitaDynamOutput(int tIndex, DynamMat& dMat, Vector& ext_f, Vector *aeroForce, 
           SysState<Vector> &state, int sliceRank, double time)
{
  startTimerMemory(times->output, times->memoryOutput);

  this->fillBcxVcx(tIndex);

  // PJSA 4-9-08 ext_f passed here may not be for the correct time
  domain->pitaDynamOutput(tIndex, bcx, dMat, ext_f, *aeroForce, state.getDisp(), state.getVeloc(),
                          state.getAccel(), state.getPrevVeloc(), vcx,
                          sliceRank, time);

  stopTimerMemory(times->output, times->memoryOutput);
}

// Update bcx for time dependent prescribed displacements and velocities (previously done in computeExternalForce2)
void
SDDynamPostProcessor::fillBcxVcx(int tIndex) {
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double t = double(tIndex)*domain->solInfo().dt;
    double *userDefineDisp = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
    double *userDefineVel  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
    //cerr << "getting usdd at time " << t << " for dynamOutput\n";
    userSupFunc->usd_disp(t,userDefineDisp,userDefineVel);
    DofSetArray *dsa = domain->getDSA();
    for(int i = 0; i < claw->numUserDisp; ++i) {
      int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
      if(dof >= 0) {
        bcx[dof] = userDefineDisp[i];
        vcx[dof] = userDefineVel[i];
      }
    }
  }
} 
