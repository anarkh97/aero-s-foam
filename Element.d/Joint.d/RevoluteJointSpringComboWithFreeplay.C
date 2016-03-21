#ifdef USE_EIGEN3
#include <Element.d/Joint.d/RevoluteJointSpringComboWithFreeplay.h>
#include <Element.d/Joint.d/RevoluteJoint.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>

RevoluteJointSpringComboWithFreeplay::RevoluteJointSpringComboWithFreeplay(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new RevoluteJoint(nnloc);
  subElems[1] = new NonlinearTorsionalSpring(nnloc, 2, 1, 1, 1);
  subElems[2] = new NonlinearTorsionalSpring(nnloc, 2, 1, 1, 2);
}

int 
RevoluteJointSpringComboWithFreeplay::getTopNumber() 
{ 
  return 106; 
}
#endif
