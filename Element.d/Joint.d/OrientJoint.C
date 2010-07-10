#include <Element.d/Joint.d/OrientJoint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>

OrientJoint::OrientJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ParallelAxesConstraintType1(nnloc);
  subElems[1] = new DotConstraintType1(nnloc, 2, 1);
}

int 
OrientJoint::getTopNumber() 
{ 
  return 106; 
}
