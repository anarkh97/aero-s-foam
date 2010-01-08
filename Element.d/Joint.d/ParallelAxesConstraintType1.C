#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>

ParallelAxesConstraintType1::ParallelAxesConstraintType1(int* _nn)
{
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new DotConstraintType1(nnloc, 1, 0);
  subElems[1] = new DotConstraintType1(nnloc, 2, 0);
  initialize(2, _nn);
}

int 
ParallelAxesConstraintType1::getTopNumber() 
{ 
  return 106; 
}

