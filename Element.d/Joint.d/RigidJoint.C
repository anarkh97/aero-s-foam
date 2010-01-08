#include <Element.d/Joint.d/RigidJoint.h>
#include <Element.d/Joint.d/ConstantDistanceConstraint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

RigidJoint::RigidJoint(int* _nn)
{
  nSubElems = 4;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ConstantDistanceConstraint(nnloc);
  subElems[1] = new ParallelAxesConstraintType1(nnloc);
  subElems[2] = new DotConstraintType1(nnloc, 2, 1);
  subElems[3] = new ParallelAxesConstraintType2(nnloc);
  initialize(2, _nn);
}

int 
RigidJoint::getTopNumber() 
{ 
  return 106; 
}
