#include <Element.d/Joint.d/RevoluteJoint.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>

RevoluteJoint::RevoluteJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(nnloc);
  subElems[1] = new ParallelAxesConstraintType1(nnloc);
}

int 
RevoluteJoint::getTopNumber() 
{ 
  return 106; 
}

