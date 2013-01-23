#ifndef _CONSTRAINTFUNCTIONELEMENT_H_
#define _CONSTRAINTFUNCTIONELEMENT_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DofSet;
class GeomState;

template<template <typename S> class ConstraintFunctionTemplate>
class ConstraintFunctionElement : public MpcElement
{
  protected:
    int rotdescr; // 0: total lagrangian, 1: updated lagrangian, 2: eulerian (default)

  public:
    ConstraintFunctionElement(int, DofSet, int*, int, int=2);
    ConstraintFunctionElement(int, DofSet*, int*, int, int=2);

    void buildFrame(CoordSet&);
    void update(GeomState&, CoordSet&, double);
    void getHessian(GeomState&, CoordSet&, FullSquareMatrix&, double t);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0.0);
    double getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    double getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t);

  protected:
    virtual void getConstants(CoordSet&,
                              Eigen::Array<typename ConstraintFunctionTemplate<double>::ScalarConstantType,
                                           ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                              Eigen::Array<int,
                                           ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                              GeomState* = NULL) {}

   virtual void getInputs(Eigen::Matrix<double,ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q, 
                          CoordSet& c0, GeomState *c1 = NULL, GeomState *refState = NULL);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/MpcElement.d/ConstraintFunctionElement.C>
#endif

#endif
