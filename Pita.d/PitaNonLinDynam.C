#include <Pita.d/PitaNonLinDynam.h>
#include <Pita.d/DynamStateSet.h>
#include <Control.d/ControlInterface.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/Domain.h>
#include <Timers.d/StaticTimers.h>
#include <iostream>

extern Communicator* structCom;

PitaNonLinDynamic::PitaNonLinDynamic(Domain *d) :
  NonLinDynamic(d),
  pitaTimers("NonLinear Pita"),
  defaultPostProcessor_(*this)
{ 
  this->preProcess();
  this->computeTimeInfo();
 
  kiter = d->solInfo().kiter;
  Jratio = d->solInfo().Jratio;
  numTSonCPU = d->solInfo().numTSperCycleperCPU;

  coarseDt = this->getDt() * Jratio;
  coarseDelta = this->getDelta() * Jratio;

  numTS = int( ceil( ( this->getTotalTime() / this->getDt() ) / Jratio) );

  totalTime = numTS * coarseDt;

  baseImprovementMethod = d->solInfo().baseImprovementMethodForPita;
}

int PitaNonLinDynamic::getInitState(DynamState<double> & ds)
{
  // Dummy vectors: We do not need that information for PITA
  GenVector<double> dummy_acc(this->solVecInfo(), 0.0);
  GenVector<double> dummy_vp(this->solVecInfo(), 0.0);
  return NonLinDynamic::getInitState(ds.disp(), ds.vel(), dummy_acc, dummy_vp);
}

int PitaNonLinDynamic::getInitSeed(DynamState<double> & ds, int sliceRank)
{
  if (sliceRank <= 0)
  {
    return getInitState(ds);
  }
  else
  {
    domain->initDispVelocOnTimeSlice(ds.disp(), ds.vel(), sliceRank);
    if (userSupFunc) // Get disp/velo when a user supplied control function has been specified
    {
      double sliceTime = domain->solInfo().getTimeStep() * domain->solInfo().Jratio * sliceRank;
      double *ctrdisp = (double *) alloca(sizeof(double)*claw->numSensor);
      double *ctrvel  = (double *) alloca(sizeof(double)*claw->numSensor);
      double *ctracc  = (double *) alloca(sizeof(double)*claw->numSensor);
      GenVector<double> dummy_acc(this->solVecInfo(), 0.0);
      extractControlData(ds.disp(), ds.vel(), dummy_acc, ctrdisp, ctrvel, ctracc);
      userSupFunc->usd_disp(sliceTime, ctrdisp, ctrvel);
    }  
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

  Connectivity *allDofs = solver->getAllDofs();
  for( iele = 0; iele < domain->numElements(); ++iele)
  {
    //int dim = kelArray[iele].dim();
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
  if (this->res != (FILE*) 0)
    fclose(this->res);
                                                                                                                                                   
  int myCPU  = structCom->myID();
  
  std::stringstream s;
  s << "residuals." << myCPU;
  this->res = fopen(s.str().c_str(), "wt");
                                                                                                                                                   
  if (this->res == (FILE *) 0)
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

void 
PitaNonLinDynamic::printNLPitaTimerFile(int CPUid)
{
  std::string fileNameString(geoSource->getCheckFileInfo()->checkfile);
  std::stringstream s;
  s << ".pitaTiming." << CPUid;
  fileNameString.append(s.str());
  std::ofstream out(fileNameString.c_str());
  if (out.fail())
  {
    fprintf(stderr, "Failed to open %s\n", fileNameString.c_str());
    return;
  }
  out.precision(5);
  out << std::fixed << pitaTimers;
  out.close();
}

// class PitaNLDynamOutput
PitaNonLinDynamic::PitaPostProcessor::PitaPostProcessor(PitaNonLinDynamic & probDesc) :
  probDesc_(probDesc),
  sliceRank_(-1)
{
}
                                                                                                                                                                                                     
PitaNonLinDynamic::PitaPostProcessor::~PitaPostProcessor()
{
  probDesc_.closeOutputFiles();
}
                                                                                                                                                                                                     
void
PitaNonLinDynamic::PitaPostProcessor::sliceRank(int rank)
{
  probDesc_.closeOutputFiles();
  if (rank >= 0)
    probDesc_.openOutputFiles(rank);
  sliceRank_ = rank;
}
