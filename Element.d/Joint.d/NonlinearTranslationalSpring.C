#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

NonlinearTranslationalSpring::NonlinearTranslationalSpring(int* _nn, int _axis, int _node)
 : DotType2ConstraintElement(_nn, _axis, _node)
{}

void
NonlinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p);
  myProp = true;
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
}

void 
NonlinearTranslationalSpring::buildFrame(CoordSet& cs)
{
  DotType2ConstraintElement::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring::update(GeomState& gState, CoordSet& cs, double t)
{
  DotType2ConstraintElement::update(gState, cs, t);
  rhs.r_value += sp0;
}
#endif
