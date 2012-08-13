#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Joint.d/exp-map.h>

NonlinearTranslationalSpring::NonlinearTranslationalSpring(int* _nn, int _axis)
 : DotConstraintType2(_nn, _axis)
{}

void
NonlinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  prop = (_myProp) ? p : new StructProp(*p);
  myProp = true;
  if(axis == 0) prop->penalty = prop->kx;
  else if(axis == 1) prop->penalty = prop->ky;
  else prop->penalty = prop->kz;
  prop->lagrangeMult = false;
}

void 
NonlinearTranslationalSpring::buildFrame(CoordSet& cs)
{
  DotConstraintType2::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring::update(GeomState& gState, CoordSet& cs, double t)
{
  DotConstraintType2::update(gState, cs, t);
  rhs.r_value += sp0;
}
