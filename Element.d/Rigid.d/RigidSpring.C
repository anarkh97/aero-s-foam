#include <Element.d/Rigid.d/RigidSpring.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/OrientJoint.h>

RigidSpring::RigidSpring(int* _nn)
{
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(indices);
  subElems[1] = new OrientJoint(indices);
  initialize(2, _nn);
}

int
RigidSpring::getTopNumber()
{
  return 102;
}

