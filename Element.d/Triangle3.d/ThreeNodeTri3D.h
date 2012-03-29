#ifndef _THREENODETRI3D_H_
#define _THREENODETRI3D_H_

#include <Element.d/Element.h>

class ThreeNodeTri3D : virtual public Element
{
 protected:
  int nn[3];

 public:
  ThreeNodeTri3D(int*);

  void renum(int *);

  int  numNodes();
  int* nodes(int * = 0);

  int  numDofs();
  int* dofs(DofSetArray &, int *p=0);
  void markDofs(DofSetArray &);

  FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
  FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
  double getMass(CoordSet&);
  void getGravityForce(CoordSet&, double *gravity, Vector& f, int gravflg, GeomState *gs);

  void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                   int surface=0, double *ndTemps=0, double ylayer=0, double zlayer=0, int avgnum=0);
  void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                    int surface=0, double *ndTemps=0);

  void computeDisp(CoordSet&cs, State &state, const InterpPoint &, double*res, GeomState *gs);
  void getFlLoad(CoordSet &, const InterpPoint &,  double *flF, double *resF, GeomState *gs=0);

  void computePressureForce(CoordSet& cs, Vector& elPressureForce, 
                            GeomState *geomState = 0, int cflg = 0, double t = 0);

  int getTopNumber();

  PrioInfo examine(int sub, MultiFront *);

  Corotator *getCorotator(CoordSet &, double *, int , int);
};

#endif
