#ifndef _RIGIDMPCSPRING_H_
#define _RIGIDMPCSPRING_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

class RigidMpcSpring : public Element, public Corotator
{
  int nn[2];
  LMPCons* mpc[6];
  double lambda[6];
  int glMpcNb[6]; 

 public:
  RigidMpcSpring(int *);
  ~RigidMpcSpring() { /* nothing to delete */ };

  void computeMPCs(CoordSet &cs, int &lmpcnum);
  void updateMPCs(GeomState &gState);
  bool isRigidMpcElement() { return true; }
  void setMpcForces(double *mpcForces) { for(int i=0; i<6; ++i) lambda[i] = mpcForces[glMpcNb[i]]; }

  Corotator *getCorotator(CoordSet &, double *, int, int) { return this; }
  void setProp(StructProp *p) { };

  void renum(int *);

  FullSquareMatrix stiffness(CoordSet &, double *, int = 1);
  FullSquareMatrix massMatrix(CoordSet &, double *, int = 1);

  void markDofs(DofSetArray &);
  int* dofs(DofSetArray &, int * = 0);
  int numDofs();
  int numNodes();
  int* nodes(int * = 0);

  void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *);
  int  getTopNumber() { return 101; }
  bool isSafe() { return true; }

};
#endif
