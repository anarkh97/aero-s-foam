#ifndef _NL_DYNAMPROBTYPE_H_
#define _NL_DYNAMPROBTYPE_H_

#include <stdio.h>
#include <stdlib.h>
#include <Driver.d/StateUpdater.h>

template < class OpSolver, 
           class VecType, 
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate = IncrUpdater<ProblemDescriptor,VecType, GeomType> >
class NLDynamSolver {
     ProblemDescriptor *probDesc;
     typename StateUpdate::RefState *refState;
     typename StateUpdate::StateIncr *stateIncr;

     FILE *vel;
     double beta, gamma, alphaf, alpham;
   public:

     // Constructor
     NLDynamSolver(ProblemDescriptor *PrbD) 
       { probDesc = PrbD; }

     void solve();
};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/NLDynamProbType.C>
#endif

#endif
