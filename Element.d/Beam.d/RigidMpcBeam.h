#ifndef _RIGIDMPCBEAM_H_
#define _RIGIDMPCBEAM_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

class RigidMpcBeam : public Element, public Corotator
{
  int nn[2];
  double len0;
  LMPCons* mpc[6];
  double lambda[6];
  int glMpcNb[6]; 

 public:
  RigidMpcBeam(int *);
  ~RigidMpcBeam() { /* nothing to delete */ };

  virtual void computeMPCs(CoordSet &cs, int &lmpcnum);
  virtual void updateMPCs(GeomState &gState);
  bool isRigidMpcElement(const DofSet &ds = DofSet::nullDofset, bool forAllNodes=false)
  { return  ds == DofSet::nullDofset || ds == (DofSet::XYZdisp | DofSet::XYZrot); }
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
