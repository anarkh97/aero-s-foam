#ifndef _VISCDAMP3D_H_
#define _VISCDAMP3D_H_
#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp.h>
#include <Math.d/FullSquareMatrix.h>

template<int nNodes>
class ViscDamp3D: public ViscDamp
{
 private:
  const static int    nDim;
  int nn[nNodes];

 public:
  ViscDamp3D(int*);
  Element* clone();

  // accessors
  int numDofs()  { return nDim*nNodes; }
  int numNodes() { return nNodes; }
  int dim()      { return nDim; }
  int* dofs(DofSetArray &, int *p=0);
  int* nodes(int * = 0);
  double* getMidPoint(CoordSet &cs);
  FullSquareMatrix dampingMatrix(CoordSet&, double* cel, int cmflg=1);
  int getTopNumber() { return 2; }

  // actions
  void markDofs(DofSetArray &);
  void renum(int*);

  int    chkOptInf(CoordSet&);
  double getGradMass(CoordSet&,CoordSet&) { return 0.0; }
  void   gradstiffness(CoordSet&,CoordSet&,FullSquareMatrix&,int = 0);
  void   gradMassMatrix(CoordSet&,CoordSet&,FullSquareMatrix&, double mratio=1);
  void   gradDampMatrix(CoordSet&,CoordSet&,FullSquareMatrix&, double freq=0.0);
  // decomposer
  bool isStart() { return false; }  
  PrioInfo examine(int sub, MultiFront *mf);

 private:
  Vector avg_normal(CoordSet&);
  void   orth_basis(CoordSet&, Vector&, Vector&, Vector&);
  void   transpose(Vector&, Vector&, Vector&, Vector&, Vector&, Vector&);
  Vector davg_normal(CoordSet&, CoordSet&);
 protected:
  double area(CoordSet&);
  double darea(CoordSet&, CoordSet&);
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp3D.C>
#endif

#endif
#endif
