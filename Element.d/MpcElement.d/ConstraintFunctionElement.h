#ifndef _CONSTRAINTFUNCTIONELEMENT_H_
#define _CONSTRAINTFUNCTIONELEMENT_H_

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/Function.d/SpaceDerivatives.h>

class DofSet;
class GeomState;

template<template <typename S> class ConstraintFunctionTemplate>
class ConstraintFunctionElement : public MpcElement
{
  protected:
#if (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
    template <typename S>
      using ConstraintJacobian = Simo::Jacobian<S,ConstraintFunctionTemplate>;
    template <typename S>
      using ConstraintJacobianVectorProduct = Simo::JacobianVectorProduct<S,ConstraintFunctionTemplate>;
#endif
  public:
    ConstraintFunctionElement(int, DofSet, int*, int);
    ConstraintFunctionElement(int, DofSet*, int*, int);

    void buildFrame(CoordSet&);
    void update(GeomState*, GeomState&, CoordSet&, double);
    void getHessian(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0.0);
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double);
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&s, double);

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
