#include <Element.d/Joint.d/RevoluteActuator.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1a.h>

RevoluteActuator::RevoluteActuator(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(nnloc);
  subElems[1] = new ParallelAxesConstraintType1(nnloc);
  subElems[2] = new DotConstraintType1a(nnloc, 2, 1);
}

int 
RevoluteActuator::getTopNumber() 
{ 
  return 106; 
}

