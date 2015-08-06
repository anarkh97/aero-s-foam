#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

NonlinearTranslationalSpring::NonlinearTranslationalSpring(int* _nn, int _axis)
 : DotType2ConstraintElement(_nn, _axis)
{}

void
NonlinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p);
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
  DotType2ConstraintElement::setProp(prop, true);
}

void 
NonlinearTranslationalSpring::buildFrame(CoordSet& cs)
{
  DotType2ConstraintElement::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  DotType2ConstraintElement::update(refState, gState, cs, t);
  rhs.r_value += sp0;
}
#endif
