#include <Element.d/Joint.d/PrismaticJointSpringCombo.h>
#include <Element.d/Joint.d/PrismaticJoint.h>
#include <Element.d/Joint.d/TranslationalSpring.h>

PrismaticJointSpringCombo::PrismaticJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new PrismaticJoint(nnloc);
  subElems[1] = new TranslationalSpring(nnloc);
}

int 
PrismaticJointSpringCombo::getTopNumber() 
{ 
  return 106; 
}
