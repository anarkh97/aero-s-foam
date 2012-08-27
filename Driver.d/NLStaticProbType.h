#ifndef _NL_STATICPROBTYPE_H_
#define _NL_STATICPROBTYPE_H_

#include<Driver.d/StateUpdater.h>

template < class OpSolver, 
           class VecType, 
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate = IncrUpdater<ProblemDescriptor, VecType, GeomType> >
class NLStaticSolver {
     ProblemDescriptor *probDesc;
     typename StateUpdate::RefState *refState;
     typename StateUpdate::StateIncr *stateIncr;
   public:

     // Constructor
     NLStaticSolver(ProblemDescriptor *PrbD) 
       { probDesc = PrbD; }

     void solve();
     void arclength();
     int  newton( VecType& force, VecType& residual, VecType &glResid,
          VecType& elementInternalForce, 
          OpSolver* solver, typename StateUpdate::RefState *refState,
          GeomType* geomState, int& numIter, double lambda=1.0, int step = 1);

     void extendedNewton(GeomType &u, GeomType &u0, VecType &dU, double &lambda, 
                double deltaLambda,
                double &deltaS, double w, int &numExtIter, OpSolver* solver,
                VecType& force, VecType& residual, VecType& glRes,
                VecType& arcLenResid,
                double forceNorm, VecType& elementInternalForce, VecType& pVec, int step = 1);

    // HB
    void predictorStep(GeomType &u, GeomType &u0, VecType &dU, double &lambda, double &deltaLambda,
                       double &deltaS, double &deltaS0, double w, OpSolver* solver,
                       VecType& force, VecType& residual, VecType &totRes, VecType& elementInternalForce,
                       VecType& duds, int step=1);

};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/NLStaticProbType.C>
#endif

#endif
