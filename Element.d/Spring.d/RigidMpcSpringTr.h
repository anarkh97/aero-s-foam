#ifndef _RIGIDMPCSPRINGTR_H_
#define _RIGIDMPCSPRINGTR_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

class RigidMpcSpringTr : public Element, public Corotator
{
  int nn[2];
  double len0;
  LMPCons **mpc;
  double *lambda;
  int *glMpcNb; 
  int nmpc;

 public:
  RigidMpcSpringTr(int *);
  ~RigidMpcSpringTr() { if(nmpc) { delete [] mpc; delete [] lambda; delete [] glMpcNb; }  };

  void setProp(StructProp *p, bool _myProp = true);

  void computeMPCs(CoordSet &cs, int &lmpcnum);
  void updateMPCs(GeomState &gState);
  bool isRigidMpcElement(const DofSet &ds = DofSet::nullDofset, bool forAllNodes=false)
     { return   ds == DofSet::nullDofset || ds == DofSet::XYZdisp; }
  void setMpcForces(double *mpcForces) { for(int i=0; i<nmpc; ++i) lambda[i] = mpcForces[glMpcNb[i]]; }

  Corotator *getCorotator(CoordSet &, double *, int, int) { return this; }

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
