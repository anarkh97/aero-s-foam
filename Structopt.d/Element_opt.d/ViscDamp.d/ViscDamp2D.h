#ifndef _VISCDAMP2D_H_
#define _VISCDAMP2D_H_
#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp.h>
#include <Math.d/FullSquareMatrix.h>

#define ViscDamp2D_nNodes 2

class ViscDamp2D: public ViscDamp
{
 private:
  const static int    nNodes;
  const static int    nDim;

  int nn[ViscDamp2D_nNodes];

 public:
  ViscDamp2D(int*);
  Element* clone();

  // accessors
  int numDofs()  { return nDim*nNodes; }
  int numNodes() { return nNodes; }
  int dim()      { return nDim; }
  int* dofs(DofSetArray &, int *p=0);
  int* nodes(int * = 0);
  double* getMidPoint(CoordSet &cs);
  FullSquareMatrix dampingMatrix(CoordSet&, double* cel,int cmflg=1);
  int getTopNumber() { return 1; }

  // actions
  void markDofs(DofSetArray &);
  void renum(int*);

  int    chkOptInf(CoordSet&);
  double getGradMass(CoordSet&,CoordSet&) { return 0.0; }
  void   gradstiffness(CoordSet&,CoordSet&,FullSquareMatrix&,int = 0);
  void   gradMassMatrix(CoordSet&,CoordSet&,FullSquareMatrix&, double mratio=1);
  void   gradDampMatrix(CoordSet&,CoordSet&,FullSquareMatrix&, double freq=0.0);

  // dec functions
  PrioInfo examine(int sub, MultiFront *mf);
  bool isStart() { return false; }  

 protected:
  double area(CoordSet&);
  double darea(CoordSet&, CoordSet&);
};
#endif
#endif
