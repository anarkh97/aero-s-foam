#include <Element.d/Joint.d/CylindricalJointSpringCombo.h>
#include <Element.d/Joint.d/CylindricalJoint.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

CylindricalJointSpringCombo::CylindricalJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new CylindricalJoint(nnloc);
  subElems[1] = new NonlinearTorsionalSpring(nnloc, 2, 1);
  subElems[2] = new NonlinearTranslationalSpring(nnloc, 0);
}

int 
CylindricalJointSpringCombo::getTopNumber() 
{ 
  return 106; 
}

