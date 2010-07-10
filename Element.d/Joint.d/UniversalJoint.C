#include <Element.d/Joint.d/UniversalJoint.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/DotConstraintType1.h>

UniversalJoint::UniversalJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(nnloc);
  subElems[1] = new DotConstraintType1(nnloc, 2, 1);
}

int 
UniversalJoint::getTopNumber() 
{ 
  return 106; 
}

