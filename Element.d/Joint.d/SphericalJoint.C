#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/LinearConstraintType1.h>

SphericalJoint::SphericalJoint(int* _nn)
{
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  DofSet xx[2] = { DofSet::Xdisp,  DofSet::Xdisp };
  DofSet yy[2] = { DofSet::Ydisp,  DofSet::Ydisp };
  DofSet zz[2] = { DofSet::Zdisp,  DofSet::Zdisp };
  subElems[0] = new LinearConstraintType1(nnloc, xx);
  subElems[1] = new LinearConstraintType1(nnloc, yy);
  subElems[2] = new LinearConstraintType1(nnloc, zz);
  initialize(2, _nn);
}

int 
SphericalJoint::getTopNumber() 
{ 
  return 106; 
}

