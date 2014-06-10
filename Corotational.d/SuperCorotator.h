#ifndef _SUPERCOROTATOR_H_
#define _SUPERCOROTATOR_H_

#include <Corotational.d/Corotator.h>
#include <Element.d/SuperElement.h>

class SuperCorotator : public Corotator
{
  int nSubElems;
  Corotator **subElemCorotators;
  SuperElement *superElem;
  FullM *origK;
  double **sub_vld;
  double **sub_vlr;

 public:
  SuperCorotator(SuperElement *_superElem); 
  virtual ~SuperCorotator();
  
  void setSubCorotator(int i, Corotator *subCorotator)
     { subElemCorotators[i] = subCorotator; }
  double* getPreviouslyExtractedSubDeformations(int i) { return (sub_vld) ? sub_vld[i] : 0; }
  double* getPreviouslyExtractedSubRigidBodyMotion(int i) { return (sub_vlr) ? sub_vlr[i] : 0; }

  void getStiffAndForce(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t);
  void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t);
  void getDExternalForceDu(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f);
  void getInternalForce(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t);
  void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t);
  void getExternalForce(GeomState &geomState, CoordSet &cs,  double *f);

  void formGeometricStiffness(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f);
  double* getOriginalStiffness();
  void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag);
  void getNLVonMises(Vector &stress, Vector &weight, GeomState &geomState, CoordSet &cs, int strInd);
  void getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState, CoordSet &cs, int strInd);
  double getElementEnergy(GeomState &geomState, CoordSet &cs);
  double getDissipatedEnergy(GeomState &geomState, CoordSet &cs);
  void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs, double *vlr);
  void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);

  void getResidualCorrection(GeomState &gs, double *r);
  void initMultipliers(GeomState& c1);
  void updateMultipliers(GeomState& c1);
  double getError();
};

#endif
