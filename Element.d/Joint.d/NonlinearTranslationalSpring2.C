#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTranslationalSpring2.h>

NonlinearTranslationalSpring2::NonlinearTranslationalSpring2(int* _nn, int _axis)
 : DotType2v2ConstraintElement(_nn, _axis)
{}

void
NonlinearTranslationalSpring2::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p);
  myProp = true;
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
}

void 
NonlinearTranslationalSpring2::buildFrame(CoordSet& cs)
{
  DotType2v2ConstraintElement::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring2::update(GeomState& gState, CoordSet& cs, double t)
{
  DotType2v2ConstraintElement::update(gState, cs, t);
  rhs.r_value += sp0;
}
#endif
