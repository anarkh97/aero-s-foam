#ifndef _DISTRTIMEDECOMPSOLVER_H_
#define _DISTRTIMEDECOMPSOLVER_H_

#include <Pita.d/PitaTimers.h>
#include <Pita.d/TimeSlice.h>
#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>
#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>
//#include <cstdlib>
//#include <algorithm>


template <
  class DynOps,
  class VecType,
  class PostProcessor,
  class ProblemDescriptor,
  class InfoSize> class LinTimeSlice;

template < 
  class VecType, 
  class ProblemDescriptor> class NewmarkWorkVec;

template < 
  class VecType,
  class DynOps,
  class InfoSize> class StateSet;

class PitaDynamMat;
class PitaTimers;
class Communicator;
class Connectivity;

template < class VecType > class SysState;

//-------------------------------------------------------------------------------
template < 
    class DynOps, 
    class VecType,
    class PostProcessor,
    class ProblemDescriptor, 
    class InfoSize>
class DistrTimeDecompSolver {

  protected:   

    // for Time management
    double Tinitial, Tfinal;
    double dt, Dt;
 
    int InitTimeIndex;
    int nTS;              // number of TimeSlice
    int numTS;            // number of TS in the CPU
    int Jratio;           // ratio coarse/fine time-grids
    int kiter;            // max number of time parallel iterations

    PitaTimers *Ptimes;

    // space data 
    InfoSize info;
    int probsize;           // total length of the problem

    // for CPU management  
    int numCPU, myCPU;
    int numTSperCycleperCPU;
    int numTotalActiveTS;
    int locidTS;
    int x;                  //x=1+numTSperCycleperCPU

    int *cycleTSid;
    int *infoCPU;

    bool *TSflag;
    bool active;
    bool nvStep0;
    bool NoForce; 
    bool ConstForce;
    bool CkCoarse;

    Connectivity *TStoCPU, *CPUtoTS;

    // for exchange
    Communicator *timeCom; 

    double *ptBuffeur;
    double *buf;

    int *numdata,*position;
    int *numdata_Info, *position_Info;

    // for computation
    double atol;
    double rtol;
 
    double normDisp_ref, normVeloc_ref;

    // TimeSlices 
    LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *arrayTimeSlices; // array of TS
    LinTimeSlice<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *ptCoarseGrid;    // pointer on TS=[Ti , ... ,Tf]

    // for solve
    ProblemDescriptor *probDesc;
    PostProcessor *postProcessor;

    SysState<VecType> *ptcurState;    // Initial condition

    NewmarkWorkVec<VecType,ProblemDescriptor> *ptworkVec;

    VecType *ptaeroForce;
    VecType *ptconstForce;

    DynOps *ptdynOps;

    // Storage array               
    StateSet<VecType,DynOps,InfoSize> *ptSeedState;   // Yi=yi(Ti)         
    StateSet<VecType,DynOps,InfoSize> *ptPropState;   // yi(Ti+1)
    StateSet<VecType,DynOps,InfoSize> *ptJumpState;   // yi-1(Ti)-Yi  
    StateSet<VecType,DynOps,InfoSize> *ptBfine;       // Bi
    StateSet<VecType,DynOps,InfoSize> *ptStep0;       // Yi for iter=0

    StateSet<VecType,DynOps,InfoSize> *ptProptmp;     // yi(Ti+1) before exchange 
    StateSet<VecType,DynOps,InfoSize> *ptCk;          // Cki
    StateSet<VecType,DynOps,InfoSize> *ptSk;          // Base (Yi)i
    StateSet<VecType,DynOps,InfoSize> *ptAmplSk;      // Amplification of this base
 
    StateSet<VecType,DynOps,InfoSize> *ptKd_Mv;       // Store (K Si_d) and (M Si_v)

    VecType *ptZero;

    int aeroAlg;
 
  public:

    DistrTimeDecompSolver(ProblemDescriptor *_probDesc);
    DistrTimeDecompSolver(){}
    virtual ~DistrTimeDecompSolver();

    // accessors
    int  getlocRank(int SliceRank);                                 //return the local rank of a TS in its CPU 
    void get_Seed(VecType &d, VecType &v, int x);

    SysState<VecType> *getptCurState()  { return ptcurState;  }   
    PitaTimers &getTimers()             { return *Ptimes;     }

    // others methods
    void buildTSandConnect();
    void getConnect();
    void printConnect(Connectivity *ptc);
    void printInfoCPU();

    void exchange_PropValue();
    void exchange_BiValue();
    void exchange_InfoCPU();
    void Ck_computation(int cycle);

    void getnumCPU();
    void arrayInitialization();
    void stopIfJumpSmall();
    void jumpEvaluation();
    void normRefcomputation();
    void newInitCond(int cycle);
    void getActiveTSidAndBuffeur();
    void myNextActiveTS();
    void vector_Initialize(bool first=false);
    void solve_PITA_linearDynam();
    void implicitNewmark_PITA(double alphaf, double alpham, double beta, double gamma, int optFlag, int *stcopt);

};


#ifdef _TEMPLATE_FIX_
#include <Pita.d/DistrTimeDecompSolverCore.C>
#include <Pita.d/DistrTimeDecompSolver.C>
#endif

#endif

