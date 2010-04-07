#ifndef _MDDISTRTIMEDECOMPSOLVER_H_
#define _MDDISTRTIMEDECOMPSOLVER_H_

class Communicator;

template <
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
class DistrTimeDecompSolver; 

//-------------------------------------------------------------------------------
template < 
    class DynOps, 
    class VecType, 
    class PostProcessor, 
    class ProblemDescriptor,
    class InfoSize>
class MDDistrTimeDecompSolver : public DistrTimeDecompSolver <DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> {

  private:   
  
   int numSpaceMPIProc;

   Communicator *timeCom;   // to communicate with MPI Process that have the same TS
 
  public:

   MDDistrTimeDecompSolver(ProblemDescriptor *_probDesc);
   MDDistrTimeDecompSolver(){}
   ~MDDistrTimeDecompSolver();


};


#ifdef _TEMPLATE_FIX_
#include <Pita.d/MDDistrTimeDecompSolver.C>
#endif

#endif

