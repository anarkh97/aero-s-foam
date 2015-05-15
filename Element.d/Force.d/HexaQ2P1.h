#ifndef _HEXAQ2P1_H_
#define _HEXAQ2P1_H_

#include <Element.d/Function.d/StrainEnergy.d/ThreeFieldStrainEnergyFunction.h>
#include <Element.d/Function.d/Shape.d/Hex20LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Linear.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Force.d/MixedFiniteElement.h>

#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
template <typename S>
using HexaQ2P1ThreeFieldStrainEnergyFunction 
      = Simo::ThreeFieldStrainEnergyFunction<S, Hex20LagrangePolynomialShapeFunction,
                                             LinearShapeFunction, LinearShapeFunction,
                                             RepeatedQuadratureRule<double,GaussLegendre,3,Eigen::Vector3d> >;

class HexaQ2P1 : public MixedFiniteElement<HexaQ2P1ThreeFieldStrainEnergyFunction>
{
  public:
    static const DofSet NODALDOFS[20];
    HexaQ2P1(int* _nn);

    int getTopNumber();
    PrioInfo examine(int sub, MultiFront *mf);
    int getQuadratureOrder(); 
};
#endif
#endif
