#ifndef _CONSTRAINTFUNCTIONELEMENT_H_
#define _CONSTRAINTFUNCTIONELEMENT_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DofSet;

template<template <typename S> class ConstraintFunctionTemplate>
class ConstraintFunctionElement : public MpcElement
{
  public:
    ConstraintFunctionElement(int, DofSet, int*, int);
    ConstraintFunctionElement(int, DofSet*, int*, int);

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
                                           ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst) {}

};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/MpcElement.d/ConstraintFunctionElement.C>
#endif

#endif
