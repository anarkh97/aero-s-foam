#ifndef _TWONODETRUSSRIGIDMPC_H_
#define _TWONODETRUSSRIGIDMPC_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

class TwoNodeTrussRigidMpc : public Element, public Corotator
{
  int nn[2];
  double len0;
  LMPCons *mpc;
  double lambda;
  int glMpcNb; 

 public:
  TwoNodeTrussRigidMpc(int *);
  ~TwoNodeTrussRigidMpc() { /* nothing to delete */ };

  virtual void computeMPCs(CoordSet &cs, int &lmpcnum);
  virtual void updateMPCs(GeomState &gState);
  bool isRigidMpcElement() { return true; }
  void setMpcForces(double *mpcForces) { lambda = mpcForces[glMpcNb]; }

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
  bool isSafe() { return false; }

};
#endif
