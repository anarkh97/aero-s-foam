#include <Driver.d/Domain.h>
#include <Utils.d/Connectivity.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::fprintf;

extern Communicator *structCom;
//-------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>
                       ::DistrTimeDecompSolver(ProblemDescriptor * _probDesc)
{
  probDesc = _probDesc;

  Ptimes = new PitaTimers;

  // TS management
  Domain *ptDom = probDesc->getDomain();

  dt       = ptDom->solInfo().getTimeStep();
  Tinitial = ptDom->solInfo().initialTime;
  Tfinal   = ptDom->solInfo().tmax;
  Jratio   = ptDom->solInfo().Jratio;
  kiter    = ptDom->solInfo().kiter;

  Dt = Jratio*dt;

  numTSperCycleperCPU = ptDom->solInfo().numTSperCycleperCPU;
  InitTimeIndex       = ptDom->solInfo().initialTimeIndex;

  // pointer for solve
  arrayTimeSlices = 0;  
  ptCoarseGrid    = 0;
  ptworkVec  = 0;  ptconstForce = 0;
  ptcurState = 0;  ptaeroForce  = 0;  
  ptdynOps   = 0;

  // pointer for storage
  ptSeedState = 0; ptProptmp = 0; 
  ptPropState = 0; ptStep0   = 0;
  ptJumpState = 0; ptKd_Mv   = 0; 

  ptCk = 0; ptBfine  = 0; 
  ptSk = 0; ptAmplSk = 0; 
 
  // pointer for CPU amnagement and exchange
  timeCom = 0;

  TStoCPU = 0; TSflag    = 0;
  CPUtoTS = 0; cycleTSid = 0; 

  ptBuffeur = 0; buf = 0; infoCPU = 0;

  numdata  = 0; position_Info = 0; 
  position = 0; numdata_Info  = 0; 

  ptZero = 0;

  normDisp_ref  = 0.0; atol = 1e-10;
  normVeloc_ref = 0.0; rtol = 1e-4;

  if (geoSource->getNewStep0()) nvStep0 = true;
  else nvStep0 = false; 
  
  if (ptDom->solInfo().ConstForcePita) ConstForce = true; 
  else ConstForce = false;

  if (ptDom->solInfo().NoForcePita) NoForce = true;
  else NoForce = false;

  if (ptDom->solInfo().CkCoarse) CkCoarse = true;
  else CkCoarse = false;

}

//-------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::~DistrTimeDecompSolver()
{
  if (ptworkVec)     { delete ptworkVec;    ptworkVec    = 0; }
  if (ptcurState)    { delete ptcurState;   ptcurState   = 0; }
  if (ptaeroForce)   { delete ptaeroForce;  ptaeroForce  = 0; }
  if (ptconstForce)  { delete ptconstForce; ptconstForce = 0; }

  if (ptSeedState)     { delete ptSeedState; ptSeedState = 0; }
  if (ptPropState)     { delete ptPropState; ptPropState = 0; }
  if (ptJumpState)     { delete ptJumpState; ptJumpState = 0; }

  if (ptBfine)   { delete ptBfine;   ptBfine   = 0; }
  if (ptProptmp) { delete ptProptmp; ptProptmp = 0; }
  if (ptSk)      { delete ptSk;      ptSk      = 0; }
  if (ptAmplSk)  { delete ptAmplSk;  ptAmplSk  = 0; }
  if (ptCk)      { delete ptCk;      ptCk      = 0; }
  if (ptStep0)   { delete ptStep0;   ptStep0   = 0; }
  if (ptKd_Mv)   { delete ptKd_Mv;   ptKd_Mv   = 0; }

  if (arrayTimeSlices) { delete[] arrayTimeSlices; arrayTimeSlices = 0; }
  if (ptCoarseGrid)    { delete ptCoarseGrid;      ptCoarseGrid    = 0; }

  if (ptZero)    { delete ptZero; ptZero=0; }

  if (Ptimes)    { delete Ptimes; Ptimes = 0; }

  if (TStoCPU)   { delete TStoCPU; TStoCPU = 0; }
  if (CPUtoTS)   { delete CPUtoTS; CPUtoTS = 0; }

  if (buf)       { delete[] buf;  buf = 0; }
  if (ptBuffeur) { delete[] ptBuffeur; ptBuffeur = 0; }
  if (numdata)   { delete[] numdata;   numdata   = 0; }
  if (position)  { delete[] position;  position  = 0; }

  if (infoCPU)       { delete[] infoCPU;       infoCPU       = 0; }
  if (position_Info) { delete[] position_Info; position_Info = 0; }
  if (numdata_Info)  { delete[] numdata_Info;  numdata_Info  = 0; }

  if (cycleTSid)     { delete[] cycleTSid; cycleTSid = 0; }

  if (TSflag)  { delete[] TSflag; TSflag=0; }

  if (ptdynOps && !CkCoarse && nvStep0) { delete ptdynOps; ptdynOps=0; }
  //If (!CkCoarse && nvStep0) ptdynOps==PitaDynamMat. It was build with new 
  //so it has to be deleted. 
  //If (CkCoarse || !nvStep0) ptdynOps==(PitaDynamMat*)DynamMat. DynamMat
  //was build with a new in buildOps, ptdynOps is a only a copy. I don't 
  //know if it has to be deleted.

}

//-------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>
                  ::get_Seed(VecType &d, VecType &v, int x)
{
    d = ptSeedState->getDisp(x);
    v = ptSeedState->getVel(x);
}

//-------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
int DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>
                ::getlocRank(int SliceRank)
{
  // find the rank of a TS in myCPU. We assume, that we know that this TS belongs 
  // to myCPU

  int locRank=0;
  while ( (*CPUtoTS)[myCPU][locRank]!=SliceRank && locRank<numTS-1 ){
    locRank++;
  }
  if ( (*CPUtoTS)[myCPU][locRank]!=SliceRank ){
    fprintf(stderr, " ... Problem in getlocRank ... \n");
    locRank=-1;
  }
  return locRank;
}

//--------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>
               ::printConnect(Connectivity *c)
{
   for (int i=0; i<c->csize(); i++){
     cout<<endl<<"pt "<<i<<" ->  ";
     for (int j=0; j<c->num(i); j++){
      cout<<(*c)[i][j]<<" ";
     }
   }
   cout<<endl;
}

//-------------------------------------------------------------------------------     
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
void DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>
               ::printInfoCPU()
{
  cout<<"infoCPU "<<endl;
  for (int i=0; i<numCPU; i++){
    cout<<infoCPU[x*i]<<"-> ";
    for (int j=1; j<infoCPU[x*i]+1; j++)
       cout<<infoCPU[x*i+j]<<" ";
    cout<<"  ";
  }
  cout<<endl;
}


