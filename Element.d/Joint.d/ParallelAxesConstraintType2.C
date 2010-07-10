#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>
#include <Element.d/Joint.d/DotConstraintType2.h>

ParallelAxesConstraintType2::ParallelAxesConstraintType2(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new DotConstraintType2(nnloc, 1);
  subElems[1] = new DotConstraintType2(nnloc, 2);
}

int 
ParallelAxesConstraintType2::getTopNumber() 
{ 
  return 106; 
}

