#include <Element.d/Joint.d/UniversalJointSpringCombo.h>
#include <Element.d/Joint.d/UniversalJoint.h>
#include <Element.d/Joint.d/TorsionalSpringType1.h>

UniversalJointSpringCombo::UniversalJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new UniversalJoint(nnloc);
  subElems[1] = new TorsionalSpringType1(nnloc, 1, 0);
  subElems[2] = new TorsionalSpringType1(nnloc, 2, 0);
}

int 
UniversalJointSpringCombo::getTopNumber() 
{ 
  return 106; 
}

