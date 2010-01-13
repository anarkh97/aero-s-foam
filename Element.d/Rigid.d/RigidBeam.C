#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Joint.d/ConstantDistanceConstraint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

RigidBeam::RigidBeam(int* _nn)
{
  nSubElems = 4;
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  subElems[0] = new ConstantDistanceConstraint(indices);
  subElems[1] = new ParallelAxesConstraintType1(indices);
  subElems[2] = new DotConstraintType1(indices, 2, 1);
  subElems[3] = new ParallelAxesConstraintType2(indices);
  initialize(2, _nn);
}
