#ifndef _NLSTATE_H_
#define _NLSTATE_H_
#include <Math.d/Vector.h>

class NLMatProbDesc;

class NLState {
   NLMatProbDesc *probDesc;
  public:
   Vector internalStates;
   Vector disp;
   Vector prescDisp;
   NLState(NLMatProbDesc * nlprd, int ninternStates, int nfree, int npresc);
   NLState(const NLState &nls);
   NLState &operator=(const NLState &);
   void update(NLState &ref, Vector &newDisp);
   void update (Vector &);
   void midpoint_step_update (Vector &, double, NLState &, Vector &);
   void get_inc_displacement (Vector &, NLState &, bool=true);
   void printNode (int);
   void print(){}

};

#endif
