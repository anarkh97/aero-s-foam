#ifndef _STATESET_H_
#define _STATESET_H_

#include <cstdio>
#include <cstdlib>

using std::fprintf;
using std::exit;

template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;

// JC: Dummy structure to enclose the metric
struct MetricMK
{
  const SparseMatrix *M, *K;
  MetricMK(SparseMatrix *_M = 0, SparseMatrix *K = 0) { M = _M; K = _M; }
  ~MetricMK() {}
};

//-------------------------------------------------------------------------------
template <
  class VecType,
  class DynOps,   
  class InfoSize>
class StateSet {

 /* 
    This class is used to store the state vector Yi=(vel;disp) used in the 
    formulation of the 2nd order hyperbolic linear problem. 

    Disp is stored on an array of pointer dispSet and Vel is stored on an 
    array of pointer velSet.

    Yi is the value obtained on the coarse time grid Ti=i*Dt with Dt=J*dt.

    StateSet has the same behaviour than a stack. New vectors are stored 
    at the end of the "stack". 

    numVectors     == number of vectors that are actually stored  
    maxnumVectors  == maximum number of vectors that can be stored
    currentVectors == index of current vector (legacy) 
 
    As StateSet is a kind of stack, if we want to change a vector stored before 
    the end, we have to change the cursor. To do this, we use adjustnumVectors(). 
    This method changes currentVector. 
 
    DynOps is a (dummy) container for the matrices M and K
    DynOps == PitaDynamMat for linear dynamics
    DynOps == MetricMK (defined above) for non-linear dynamics

    InfoSize is used to build the pointer for the new vectors. 
    * VecType(Infosize) with :

    InfoSize == int for Single Space Domain 
    InfoSize == DistrInfo for Multi Space Domain (not used yet)

 */

 /*  
    len == total length of a vector. this data is used to check if vectors
    stored have a correct size. This test is used in the Single Space Domain 
    case. It might be changed for Multi Space Doamain case.
 */

 private:
    
    static const int adjustMaxStep = 10; 
 
 protected:

    int numVectors; 
    int maxnumVectors;
    int currentVector; // JC: Legacy (internal index, should be replaced by an external index)

    InfoSize info;        
    int len;              // total length of vectors stored

    VecType **dispSet;    // array of pointer of Vectors
    VecType **velSet;

 public:

    StateSet();
    StateSet(const StateSet &ss);
    StateSet(int maxnum, const InfoSize &inf);
    ~StateSet();

    // Global operations
    void copySet(const StateSet &ss);
    void mergeSet(const StateSet &ss);
    void emptySet();
    void newSet(int maxnum, const InfoSize &inf);

    // Mutators
    void setMaxNumVectors(int n);
    void setState(const VecType &disp, const VecType &vel, int i);

    // Accessors
    int getSize()          const { return len; }
    int getnumVectors()    const { return numVectors; }
    int getmaxnumVectors() const { return maxnumVectors; }
    InfoSize getInfo()     const { return info; }
  
    //JC: Not used and unsafe 
    //VecType *getptdispState(int i) const { return dispSet[i]; }
    //VecType *getptvelState(int i)  const { return velSet[i]; }

    const VecType & getDisp(int i) const {
       if (i < numVectors)
         return *(dispSet[i]);
       else
       {
         fprintf(stderr," ... Error in StateSet::getDisp : out of range ... "); 
         exit(1); //return new VecType(len);
       }
    }

    const VecType & getVel(int i) const  {
       if (i < numVectors)
         return *(velSet[i]);
       else
       {
         fprintf(stderr," ... Error in StateSet::getVel : out of range ... ");
         exit(1); //return new VecType(len);
       }
    }

    //JC: Legacy (unsafe) method
    VecType & getdispState(int i) const {
       if (i < numVectors)
         return *(dispSet[i]);
       else
       {
         fprintf(stderr," ... Error in StateSet::getdispState : out of range ... ");
         exit(1); //return new VecType(len);
       }
    }
    
    VecType & getvelState(int i) const  {
       if (i < numVectors)
         return *(velSet[i]);
       else
       {
         fprintf(stderr," ... Error in StateSet::getvelState : out of range ... ");
         exit(1); //return new VecType(len);
       }
    } 


    // --- Legacy methods (for use with currentVector index)
    // "Set[numVectors]=S[j]" numVectors++
    void put_end_StateSet(const VecType &disp, const VecType &vel);      
    void put_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i); 
    void put_end_StateSet(double *array, int num, const VecType &vel);

    // "Set[i]+=S[j]"
    void add_end_StateSet(const VecType &disp, const VecType &vel);      
    void add_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i);          
    void add_end_StateSet(const VecType &disp, double *array, int num);

    // "Set[i]-=S[j]"
    void sub_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i);

    // "Set[numVectors-1]=S[j]"
    void replace_end_StateSet(const VecType &disp, const VecType &vel);
    void replace_end_StateSet(const StateSet<VecType,DynOps,InfoSize> &S, int i);
 
    void replace_StateSet(int k, const VecType &disp, const VecType &vel);
    
    void adjustnumVectors(int i=0); 
    // --- End of legacy methods

    //JC: (Legacy) method to build base
    void buildBases(StateSet<VecType,DynOps,InfoSize> &Prop, StateSet<VecType,DynOps,InfoSize> &Sk, 
        StateSet<VecType,DynOps,InfoSize> &GFSk, DynOps &dynOps, StateSet<VecType,DynOps,InfoSize> &Kd_Mv, 
        StateSet<VecType,DynOps,InfoSize> &Bi, int numActiveTS, int *cycleTSid);
 
    void buildOrtho(StateSet<VecType,DynOps,InfoSize> &final, const DynOps &dynOps);

};

//--------------------------------------------------------------------------------------------------------------

#ifdef _TEMPLATE_FIX_
#include <Pita.d/StateSet.C>
#endif

#endif


