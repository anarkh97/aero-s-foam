#if !defined(_NLHEXAHEDRAL_H_) && defined(USE_EIGEN3)
#define _NLHEXAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>
#include <Element.d/Function.d/Shape.d/Hex8LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Hex20LagrangePolynomial.h>
#include <Element.d/Function.d/Shape.d/Hex32LagrangePolynomial.h>

typedef SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8> NLHexahedral8;
typedef SolidElementTemplate<Hex20LagrangePolynomialShapeFunction,20,27> NLHexahedral20;
typedef SolidElementTemplate<Hex32LagrangePolynomialShapeFunction,32,64> NLHexahedral32;

// Nonlinear 8-node brick element with 2x2x2 gauss points
class NLHexahedral : public SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>
{
    bool linearKinematics;
  public:
    NLHexahedral(int *nd, bool isLinKin) : linearKinematics(isLinKin), SolidElementTemplate<Hex8LagrangePolynomialShapeFunction,8,8>(nd) {}
    StrainEvaluator* getStrainEvaluator();
    int getTopNumber() override { return 117; }
};

#endif
