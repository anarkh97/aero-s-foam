#ifndef _FORCEFUNCTIONELEMENT_H_
#define _FORCEFUNCTIONELEMENT_H_

#include <Element.d/Force.d/BoundaryElement.h>

class DofSet;
class GeomState;

template<template <typename S> class VectorValuedFunctionTemplate>
class ForceFunctionElement : public BoundaryElement
{
  public:
    ForceFunctionElement(int, DofSet, int*);
    ForceFunctionElement(int, DofSet*, int*);
    ForceFunctionElement(int, DofSet*, DofSet*, int*);

    FullSquareMatrix stiffness(CoordSet&, double*, int = 1);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getInternalForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0.0);
  private:
    void getJacobian(GeomState *refState, GeomState &c1, CoordSet& c0, FullM& B, double t);

  protected:
    virtual void getConstants(CoordSet&,
                              Eigen::Array<typename VectorValuedFunctionTemplate<double>::ScalarConstantType,
                                           VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                              Eigen::Array<int,
                                           VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                              GeomState* = NULL, double = 0) {}

   virtual void getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q, 
                          CoordSet& c0, GeomState *c1 = NULL, GeomState *refState = NULL);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Force.d/ForceFunctionElementImpl.h>
#endif

#endif
