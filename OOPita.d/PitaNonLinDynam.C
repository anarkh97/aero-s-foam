#include <OOPita.d/PitaNonLinDynam.h>
#include <OOPita.d/DynamStateBasis.h>
#include <Control.d/ControlInterface.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/Domain.h>
#include <Timers.d/StaticTimers.h>
#include <iostream>
#include <algorithm>

extern Communicator* structCom;

namespace Pita {

PitaNonLinDynamic::PitaNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain)
{}

void
PitaNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  computeTimeInfo();
 
  kiter = domain->solInfo().kiter;
  Jratio = domain->solInfo().Jratio;
  numTSonCPU = domain->solInfo().numTSperCycleperCPU;

  coarseDt = getDt() * Jratio;
  coarseDelta = getDelta() * Jratio;

  numTS = int( ceil( ( getTotalTime() / getDt() ) / Jratio) );

  totalTime = numTS * coarseDt;

  baseImprovementMethod = domain->solInfo().baseImprovementMethodForPita;
}

int PitaNonLinDynamic::getInitSeedCount() const {
  return geoSource->getNewStep0() ? std::max(std::min(std::min(geoSource->getNumTSPitaIDis6(), geoSource->getNumTSPitaIVel6()), numTS), 1) : 1;
}

int PitaNonLinDynamic::getInitState(DynamState & ds)
{
  // Dummy vectors: We do not need that information for PITA
  GenVector<double> dummy_acc(solVecInfo(), 0.0);
  GenVector<double> dummy_vp(solVecInfo(), 0.0);
  return NonLinDynamic::getInitState(ds.displacement(), ds.velocity(), dummy_acc, dummy_vp);
}

int PitaNonLinDynamic::getInitSeed(DynamState & ds, int sliceRank)
{
  if (sliceRank <= 0)
  {
    return getInitState(ds);
  }
  else
  {
    domain->initDispVelocOnTimeSlice(ds.displacement(), ds.velocity(), sliceRank);
    double sliceTime = domain->solInfo().getTimeStep() * domain->solInfo().Jratio * sliceRank;
    GenVector<double> dummy_acc(solVecInfo(), 0.0);
    GenVector<double> dummy_vp(solVecInfo(), 0.0);
    updateUserSuppliedFunction(ds.displacement(), ds.velocity(), dummy_acc, dummy_vp, sliceTime);
    return 0; // Default value for int aeroAlg
  }
}

// Rebuild dynamic mass matrix and stiffness matrix (fine time-grid)
void PitaNonLinDynamic::reBuildKonly()
{
  times->rebuild -= getTime();

  int iele;

  if (kuc) kuc->zeroAll();
  if (K) K->zeroAll();

  Connectivity *allDofs = domain->getAllDOFs();
  for( iele = 0; iele < domain->numElements(); ++iele)
  {
    if (kuc) kuc->add( kelArray[iele], (*allDofs)[iele] );
    if (K) K->add( kelArray[iele], (*allDofs)[iele] );
  }

  times->rebuild += getTime();
}

// Set rotational displacements equal to zero.
void PitaNonLinDynamic::zeroRotDofs(Vector & vec) const
{
  ConstrainedDSA & cdsa = *(domain->getCDSA());
  int numNodes = cdsa.numNodes(); 
  int dofPos;
  for (int inode = 0; inode < numNodes; ++inode)
  {
      dofPos = cdsa.locate(inode, DofSet::Xrot);
      if (dofPos >= 0)
        vec[dofPos] = 0.0;
      dofPos = cdsa.locate(inode, DofSet::Yrot);
      if (dofPos >= 0)
        vec[dofPos] = 0.0;
      dofPos = cdsa.locate(inode, DofSet::Zrot);
      if (dofPos >= 0)
        vec[dofPos] = 0.0;
  }
}

void PitaNonLinDynamic::buildOps(AllOps<double> & allOps, double Kcoef, double Mcoef, double Ccoef, Rbm * rigidBodyModes)
{
  allOps.K = domain->constructDBSparseMatrix<double>();
  domain->buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, rigidBodyModes);
  K = allOps.K;
}

double PitaNonLinDynamic::energyNorm(const Vector &disp, const Vector &velo)
{
  return sqrt(energyDot(disp, velo, disp, velo)); 
}

double PitaNonLinDynamic::energyDot(const Vector &disp1, const Vector &velo1, const Vector &disp2, const Vector &velo2)
{
  Vector Kdisp(disp1.size());
  Vector Mvelo(velo1.size());
  K->mult(disp1, Kdisp);
  M->mult(velo1, Mvelo);
  return (Mvelo * velo2) + (Kdisp * disp2); 
}

void PitaNonLinDynamic::openResidualFile()
{
  if (res != (FILE*) 0)
    fclose(res);
  
  int myCPU  = structCom->myID();
  
  std::stringstream s;
  s << "residuals." << myCPU;
  res = fopen(s.str().c_str(), "wt");
  
  if (res == (FILE *) 0)
    filePrint(stderr, " *** ERROR: Cannot open residual file for CPU # %d\n", myCPU);
}

// No Aero
void
PitaNonLinDynamic::pitaDynamOutput(int timeSliceRank, GeomState* geomState, Vector& velocity,
                                   Vector& vp, double time, int step, Vector& force, Vector &aeroF)
{
  times->output -= getTime();
  domain->pitaPostProcessing(timeSliceRank, geomState, force, aeroF, time, step, velocity.data(), vcx, allCorot, melArray);
  times->output += getTime();
}

void
PitaNonLinDynamic::openOutputFiles(int sliceRank)
{
  geoSource->openOutputFilesForPita(sliceRank);
  //domain->printStatistics(); // Deactivated
}

void
PitaNonLinDynamic::closeOutputFiles()
{
  geoSource->closeOutputFiles();
}

} // end namespace Pita
