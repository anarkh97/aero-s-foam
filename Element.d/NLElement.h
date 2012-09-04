#ifndef _NLELEMENT_H_
#define _NLELEMENT_H_
#include <Element.d/Element.h>
#include <cstdio>

// Declaration of the Material Non linear element
class MatNLElement : public Element { 
   public:
     MatNLElement() {}
/*
     // numStates() returns the number of states present in the state history
     // vector for this element
     virtual int numStates() { return 0; }
     virtual void initStates(double *) {}
*/
     virtual void getStiffAndForce(Node *nodes, double *disp,
                                   double *state, FullSquareMatrix &kTan,
                                   double *force) {
       fprintf(stderr, "MatNLElement::getStiffAndForce is being called on an element "
               "for which it is not defined\n");
     }
     virtual void integrate(Node *nodes, double *dispn, double *staten,
                            double *dispnp, double *statenp,
                            FullSquareMatrix &kTan,
                            double *force, double dt=0.0) {
       int nst = numStates();
       int i;
       for(i = 0; i < nst; ++i)
         statenp[i] = staten[i];
         updateStates(nodes, dispn, dispnp, statenp);
         getStiffAndForce(nodes, dispnp, statenp, kTan, force);
     }
     virtual void getInternalForce(Node *nodes, double *disp,
                                   double *state, double *force) {
       fprintf(stderr, "MatNLElement::getInternalForce is being called on an element "
               "for which it is not defined\n");
     }
     virtual void integrate(Node *nodes, double *dispn, double *staten,
                            double *dispnp, double *statenp,
                            double *force, double dt=0.0) {
       int nst = numStates();
       int i;
       for(i = 0; i < nst; ++i)
         statenp[i] = staten[i];
         updateStates(nodes, dispn, dispnp, statenp);
         getInternalForce(nodes, dispnp, statenp, force);
     }

     virtual void updateStates(Node *nodes, double *un, double *unp, double *statenp) {
       fprintf(stderr, "MatNLElement::updateStates is being called on an element "
               "for which it is not defined\n");
     }
     virtual void getStrainTens(Node *nodes, double *dispnp, double (*result)[9]) {
       fprintf(stderr, "MatNLElement::getStrainTens is being called on an element "
               "for which it is not defined\n");
     }
     virtual void getVonMisesStrain(Node *nodes, double *dispnp, double *result) {
       fprintf(stderr, "MatNLElement::getVonMisesStrain is being called on an element "
               "for which it is not defined\n");
     }
     virtual void getStressTens(Node *nodes, double *dispn, double *staten,
                                double *dispnp, double *statenp, double (*result)[9]) {
       fprintf(stderr, "MatNLElement::getStressTens is being called on an element "
               "for which it is not defined\n");
     }
     virtual void getVonMisesStress(Node *nodes, double *dispn, double *staten,
                                    double *dispnp, double *statenp, double *result) {
       fprintf(stderr, "MatNLElement::getVonMisesStress is being called on an element "
               "for which it is not defined\n");
     }
     virtual void getEquivPlasticStrain(double *statenp, double *result) {
       fprintf(stderr, "MatNLElement::getEquivPlasticStrain is being called on an element "
               "for which it is not defined\n");
     }
};

#endif
