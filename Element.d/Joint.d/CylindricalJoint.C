#include <Element.d/Joint.d/CylindricalJoint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

CylindricalJoint::CylindricalJoint(int* _nn)
{
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ParallelAxesConstraintType1(nnloc);
  subElems[1] = new ParallelAxesConstraintType2(nnloc);
  initialize(2, _nn);
}

int 
CylindricalJoint::getTopNumber() 
{ 
  return 106; 
}
