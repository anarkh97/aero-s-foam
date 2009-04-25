#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_

template <
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize> class DistrTimeDecompSolver;

template <
    class VecType,
    class DynOps,
    class InfoSize> class StateSet;

//#include <Driver.d/DynamProbType.h>

template <typename VecType> class SysState;

template <
    class VecType,
    class ProblemDescriptor> class NewmarkWorkVec;

//-------------------------------------------------------------------------------------
/*template <
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
class GenTimeSlice {

  protected:
  
    double Ti_initial;
    double Ti_final;
    double dt;
    
    int SliceRank;               // SliceRank: {0,...,nTS-1}
    int TimeIndex;
    
    int numSteps;                // number of interval = number of time iteration
    int nTS;                     // number of TimeSlice                                                                                                                                                                                
    int numITA;
    int numTSinmyCPU;
    int locRank;
                                                                                                                                                                                    
    bool convergence;
                                                                                                                                                                                    
    ProblemDescriptor *probDesc;
    PostProcessor     *postProcessor;
                                                                                                                                                                                    
    DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *ptTimeDecompSolver;

    DynOps *ptdynOps;   

    double alphaf, alpham, beta, gamma;

  public:

    GenTimeSlice();
    ~GenTimeSlice();
                                                                                                                                                                                    
    // Accessors
    int  getSliceRank()   {return SliceRank;}
    int  getnumITA()      {return numITA;}
    bool getconvergence() {return convergence;}
                                                                                                                                                                                    
    void getParam(int &numsave, double &dtsave) {numsave=numSteps; dtsave=dt;}
                                                                                                                                                                                    
    // Others method
    void setProblemParam(DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *_pt,
                                                ProblemDescriptor *_probDesc, PostProcessor *_postProcessor);
    void setTimeParam(double _Ti_initial, double _Ti_final, double _dt, int _tIndex, int _numTimeSteps);
    void setParam(int _SliceRank, int _locRank, int _numTSinmyCPU);
                                                                                                                                                                                    
    //void changeParam(bool ck, int numsave=0, double dtsave=0.0);
                                                                                                                                                                                    
    //void initialization(VecType *ptconstForce, VecType *aeroForce, double alphaf, double alpham,
    //                    double beta,double gamma, DynOps *ptdynOps);
                                                                                                                                                                                    
    //void getInitCond(bool step0, bool ck, bool bi, int myCPU, int cycleTSid_0, StateSet<VecType,DynOps,InfoSize> *tmp=0);
    
    void converge()      {convergence=true;}
    void addnumITA()     {numITA++;}
    void print();

};
*/

//-------------------------------------------------------------------------------------
template < 
    class DynOps, 
    class VecType, 
    class PostProcessor, 
    class ProblemDescriptor,
    class InfoSize>
class LinTimeSlice { 

  protected:
                                                                                                                                                                                                     
    double Ti_initial;
    double Ti_final;
    double dt;
                                                                                                                                                                                                     
    int SliceRank;               // SliceRank: {0,...,nTS-1}
    int TimeIndex;
                                                                                                                                                                                                     
    int numSteps;                // number of interval = number of time iteration
    int nTS;                     // number of TimeSlice                                                                                                                                                                                                      
    int numITA;
    int numTSinmyCPU;
    int locRank;
                                                                                                                                                                                                     
    bool convergence;
                                                                                                                                                                                                     
    ProblemDescriptor *probDesc;
    PostProcessor     *postProcessor;
                                                                                                                                                                                                     
    DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *ptTimeDecompSolver;
                                                                                                                                                                                                     
    DynOps *ptdynOps;
                                                                                                                                                                                                     
    double alphaf, alpham, beta, gamma;
                                                                                                                                                                                                     
  public:
                                                                                                                                                                                                     
    // Accessors
    int  getSliceRank()   {return SliceRank;}
    int  getnumITA()      {return numITA;}
    bool getconvergence() {return convergence;}
                                                                                                                                                                                                     
    void getParam(int &numsave, double &dtsave) {numsave=numSteps; dtsave=dt;}
                                                                                                                                                                                                     
    // Others method
    void setProblemParam(DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> *_pt,
                                                ProblemDescriptor *_probDesc, PostProcessor *_postProcessor);
    void setTimeParam(double _Ti_initial, double _Ti_final, double _dt, int _tIndex, int _numTimeSteps);
    void setParam(int _SliceRank, int _locRank, int _numTSinmyCPU);
                                                                                                                                                                                                     
    void converge()      {convergence=true;}
    void addnumITA()     {numITA++;}
    void print();

 protected:

    SysState<VecType> *ptcurState;                  // Initial conditions
    VecType *d, *v, *a, *vp; 

    NewmarkWorkVec<VecType,ProblemDescriptor> *ptworkVec;
    
    VecType *ptaeroForce;
    VecType *ptconstForce;
  
 public:

    LinTimeSlice();
    ~LinTimeSlice();

    //accessors
    NewmarkWorkVec<VecType,ProblemDescriptor> *getptWorkVec() {return ptworkVec;}
    SysState<VecType>                        *getptCurState() {return ptcurState;}

    void changeParam(bool ck, int numsave = 0, double dtsave = 0.0);

    //others method
    void initialization(VecType *ptconstForce, VecType *aeroForce, double alphaf, double alpham,
                        double beta,double gamma, DynOps *ptdynOps);

    void getInitCond(bool step0, bool ck, bool bi, int myCPU, int cycleTSid_0, StateSet<VecType,DynOps,InfoSize> *tmp=0);
    void solve_ITA_linearDynam(bool bi, bool step0, bool ck, int myCPU, StateSet<VecType,DynOps,InfoSize> *tmp, 
                               int cycleTSid_0);
};

#ifdef _TEMPLATE_FIX_
#include <Pita.d/TimeSlice.C>
#endif

#endif


