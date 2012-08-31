#ifdef USE_EIGEN3
#include <Element.d/Joint.d/UniversalJointSpringCombo.h>
#include <Element.d/Joint.d/UniversalJoint.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>

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
  subElems[1] = new NonlinearTorsionalSpring(nnloc, 2, 0);
  subElems[2] = new NonlinearTorsionalSpring(nnloc, 1, 0);
}

void
UniversalJointSpringCombo::setProp(StructProp *p, bool myProp)
{
  SuperElement::setProp(p, myProp);

  subElems[1]->getProperty()->penalty = p->k1;
  subElems[2]->getProperty()->penalty = p->k2;
}

int 
UniversalJointSpringCombo::getTopNumber() 
{ 
  return 106; 
}
#endif
