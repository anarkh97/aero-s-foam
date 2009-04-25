#ifndef _SDDISTRTIMEDECOMPSOLVER_C_
#define _SDDISTRTIMEDECOMPSOLVER_C_
                                                                                                                       
#include <Pita.d/DistrTimeDecompSolver.h>
#include <Pita.d/SDDistrTimeDecompSolver.h>

class Communicator;

extern Communicator *structCom;

//----------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
SDDistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::SDDistrTimeDecompSolver(ProblemDescriptor * 
   _probDesc):DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>(_probDesc)
{
  this->timeCom = structCom;
  
}

//----------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
SDDistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::~SDDistrTimeDecompSolver()
{
  //Don't delete timeCom. timeCom=structCom so it will be deleted by deleting structCom.
  //This will be done in main.C
}

//----------------------------------------------------------------------------------
#endif
