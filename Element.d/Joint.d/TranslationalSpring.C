#include <Element.d/Joint.d/TranslationalSpring.h>

TranslationalSpring::TranslationalSpring(int* nn)
 : ConstantDistanceConstraint(nn)
{
}

void
TranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p);
  myProp = true;
  prop->penalty = prop->kx;
  prop->lagrangeMult = false;
}
