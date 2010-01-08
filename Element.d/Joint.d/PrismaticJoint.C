#include <Element.d/Joint.d/PrismaticJoint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

PrismaticJoint::PrismaticJoint(int* _nn)
{
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ParallelAxesConstraintType1(nnloc);
  subElems[1] = new DotConstraintType1(nnloc, 2, 1);
  subElems[2] = new ParallelAxesConstraintType2(nnloc);
  initialize(2, _nn);
}

int 
PrismaticJoint::getTopNumber() 
{ 
  return 106; 
}
