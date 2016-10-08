#ifndef _MIXEDFINITEELEMENT_H_
#define _MIXEDFINITEELEMENT_H_

#include <Element.d/Force.d/BoundaryElement.h>

class DofSet;
class GeomState;
class NLMaterial;

template<template <typename S> class ScalarValuedFunctionTemplate>
class MixedFiniteElement : public BoundaryElement
{
    enum {
      N = ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,
      NumberOfNodes = ScalarValuedFunctionTemplate<double>::NumberOfNodes,
      NumberOfDimensions = ScalarValuedFunctionTemplate<double>::NumberOfDimensions
    };
    bool first_time;
    Eigen::Matrix<double,N,1> q_copy;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> CinvBt;
    Eigen::Matrix<double,Eigen::Dynamic,1> Cinvg;
    double epsilon; // penalty parameter for incompressible materials

  protected:
    int materialType;
    NLMaterial *mat;
    int nIV; // number of internal variables (excluding lagrange multipliers for incompressible materials)
    int nLM; // number of lagrange multiplier variables for incompressible materials

  public:
    MixedFiniteElement(int, DofSet, int*);
    MixedFiniteElement(int, DofSet*, int*);

    void setMaterial(NLMaterial *_mat);
    bool isSafe() { return true; }

    FullSquareMatrix stiffness(CoordSet&, double*, int = 1);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getStiffAndForce(GeomState*, GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getHessian(GeomState *refState, GeomState *c1, CoordSet& c0, Eigen::Matrix<double,
                    ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,
                    ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates>& H, double t);

    int numStates() { return nIV+nLM; }
    void initStates(double *states) { for(int i=0; i<nIV+nLM; ++i) states[i] = 0; }

    void initMultipliers(GeomState& c1);
    double getError(GeomState& c1);
    void updateMultipliers(GeomState& c1);

    virtual int getQuadratureOrder() = 0;

  private:
    void getConstants(CoordSet &cs, Eigen::Array<typename ScalarValuedFunctionTemplate<double>::ScalarConstantType,
                      ScalarValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                      Eigen::Array<int, ScalarValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst,
                      GeomState *gs = NULL);
    void getInputs(Eigen::Matrix<double,ScalarValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q, 
                   CoordSet& c0, GeomState *c1 = NULL);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Force.d/MixedFiniteElement.C>
#endif

#endif