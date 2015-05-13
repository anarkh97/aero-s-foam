#ifndef _HEXAQ2Q1_H_
#define _HEXAQ2Q1_H_

#include <Element.d/Function.d/StrainEnergy.d/ThreeFieldStrainEnergyFunction.h>
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Hex20LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Force.d/MixedFiniteElement.h>

#if (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
template <typename S>
using HexaQ2Q1ThreeFieldStrainEnergyFunction 
      = Simo::ThreeFieldStrainEnergyFunction<S, Hex20LagrangePolynomialShapeFunction,
                                             Hex8LagrangePolynomialShapeFunction, Hex8LagrangePolynomialShapeFunction,
                                             RepeatedQuadratureRule<double,GaussLegendre,3,Eigen::Vector3d> >;

class HexaQ2Q1 : public MixedFiniteElement<HexaQ2Q1ThreeFieldStrainEnergyFunction>
{
  public:
    static const DofSet NODALDOFS[20];
    HexaQ2Q1(int* _nn);

    int getTopNumber();
    PrioInfo examine(int sub, MultiFront *mf);
    int getQuadratureOrder(); 
};
#endif
#endif
