#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>
#include <Element.d/Joint.d/ConstantDistanceConstraint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>

RigidTwoNodeTruss::RigidTwoNodeTruss(int* _nn)
  : ConstantDistanceConstraint(_nn)
{
}
