#include <Element.d/Joint.d/PinInSlotJoint.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

PinInSlotJoint::PinInSlotJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new DotConstraintType1(nnloc, 2, 1);
  subElems[1] = new DotConstraintType1(nnloc, 1, 0);
  subElems[2] = new ParallelAxesConstraintType2(nnloc);
}

int 
PinInSlotJoint::getTopNumber() 
{ 
  return 106; 
}
