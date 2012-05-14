#include <Element.d/Joint.d/PlanarJoint.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType2.h>

PlanarJoint::PlanarJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new DotConstraintType1(nnloc, 2, 1);
  subElems[1] = new DotConstraintType1(nnloc, 2, 0);
  subElems[2] = new DotConstraintType2(nnloc, 2);
}

int 
PlanarJoint::getTopNumber() 
{ 
  return 106; 
}
