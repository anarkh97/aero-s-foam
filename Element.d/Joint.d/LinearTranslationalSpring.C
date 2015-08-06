#ifdef USE_EIGEN3
#include <Element.d/Joint.d/LinearTranslationalSpring.h>

LinearTranslationalSpring::LinearTranslationalSpring(int* nn)
 : ConstantDistanceConstraint(nn)
{
}

void
LinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p);
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
  ConstantDistanceConstraint::setProp(prop, true);
}
#endif
