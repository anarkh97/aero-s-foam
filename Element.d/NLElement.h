#ifndef _NLELEMENT_H_
#define _NLELEMENT_H_
#include <Element.d/Element.h>
#include <cstdio>

// Declaration of the Material Non linear element
class MatNLElement : public Element { 
   public:
     MatNLElement() {}
     // numStates() returns the number of states present in the state history
     // vector for this element
     
     virtual int numStates() { return 0; }
     virtual void initStates(double *) {}
     virtual void getStiffAndForce(Node *nodes, double *disp,
                         double *state, FullSquareMatrix &kTan,
                         double *force) {
       fprintf(stderr, "getStiffAndForce is being called on an element"
              "for which it is not defined\n");
     }
     virtual void integrate(Node *nodes, double *dispn,  double *staten,
                                         double *dispnp, double *statenp,
                                         FullSquareMatrix &kTan,
                         double *force, double dt=0.0) {
                int nst = numStates();
                int i;
                for(i = 0; i < nst; ++i)
                  statenp[i] = staten[i];
                updateStates(nodes, statenp, dispn, dispnp);
                getStiffAndForce(nodes, dispnp, statenp, kTan, force);
           }
     virtual void updateStates(Node *nodes, double *states,
                     double *un, double *unp) {
       fprintf(stderr, "updateStates is being called on an element"
              "for which it is not defined\n");
     }
    
};

#endif
