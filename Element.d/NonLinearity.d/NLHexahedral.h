#ifndef _NLHEXAHEDRAL_H_
#define _NLHEXAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>

typedef SolidElementTemplate<Hexahedron,8, 8 > NLHexahedral8;
typedef SolidElementTemplate<Hexahedron,20,27> NLHexahedral20;
typedef SolidElementTemplate<Hexahedron,32,64> NLHexahedral32;

// Nonlinear 8-node brick element with 2x2x2 gauss points
class NLHexahedral : public SolidElementTemplate<Hexahedron,8,8>
{
    bool linearKinematics;
  public:
    NLHexahedral(int *nd, bool isLinKin) : linearKinematics(isLinKin), SolidElementTemplate<Hexahedron,8,8>(nd) {}
    StrainEvaluator* getStrainEvaluator();
    PrioInfo examine(int sub, MultiFront *);
    int getTopNumber() { return 117; }
};

#endif
