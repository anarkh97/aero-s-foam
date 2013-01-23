#ifdef USE_EIGEN3
#include <Element.d/Joint.d/LinearTranslationalSpring.h>

LinearTranslationalSpring::LinearTranslationalSpring(int* nn)
 : ConstantDistanceConstraint(nn)
{
}

void
LinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p);
  myProp = true;
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
}
#endif
