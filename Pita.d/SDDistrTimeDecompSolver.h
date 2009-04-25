#ifndef _SDDISTRTIMEDECOMPSOLVER_H_
#define _SDDISTRTIMEDECOMPSOLVER_H_

#include <Pita.d/DistrTimeDecompSolver.h>

template < 
    class DynOps, 
    class VecType, 
    class PostProcessor, 
    class ProblemDescriptor,
    class InfoSize>
class SDDistrTimeDecompSolver : public DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize> {

  public:

   SDDistrTimeDecompSolver(ProblemDescriptor *_probDesc);
   SDDistrTimeDecompSolver(){}
   ~SDDistrTimeDecompSolver();


};


#ifdef _TEMPLATE_FIX_
#include <Pita.d/SDDistrTimeDecompSolver.C>
#endif

#endif

