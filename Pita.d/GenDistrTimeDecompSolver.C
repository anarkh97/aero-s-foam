#include <Pita.d/GenDistrTimeDecompSolver.h>
#include <Pita.d/DynamStateSet.h>
#include <Comm.d/Communicator>

extern Communicator structCom;

template <typename VecType>
GenDistrTimeDecompSolver<VecType>::GenDistrTimeDecompSolver(NonLinDynamic * pbdesc)
{
  probDesc = pbdesc;
  domain = probDesc->getDomain();
}

template <typename Vectype>
void GenDistrTimeDecompSolver<VecType>::initialize() {

  domain = probDesc->getDomain();
  timeCom = structCom;
 
  dt       = ptDom->solInfo().dt;
  Tinitial = ptDom->solInfo().initialTime;
  Tfinal   = ptDom->solInfo().tmax;
  Jratio   = ptDom->solInfo().Jratio;
  kiter    = ptDom->solInfo().kiter;
  Dt = Jratio * dt;

  numCPU = timeCom->numCPUs();
  myCPU  = timeCom->myID();

  InitTimeIndex       = ptDom->solInfo().initialTimeIndex;
  maxNumTSonCPU = probDesc->solInfo().numTSperCycleperCPU;

  timeSlice = new NLTimeSlice[maxNumTSonCPU]; 

}

virtual template <typename VecType>
GenDistrTimeDecompSolver<VecType>::~GenDistrTimeDecompSolver()
{
   delete[] timeSlice;                                                                                                             
}
