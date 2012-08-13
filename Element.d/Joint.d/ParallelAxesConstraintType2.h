#ifndef _PARALLELAXESCONSTRAINTTYPE2_H_
#define _PARALLELAXESCONSTRAINTTYPE2_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.3
// Constraint Building Block type BB3, also known as a straight line point follower constraint
// two constrained translational dofs

class ParallelAxesConstraintType2 : public SuperElement
{
  public:
    ParallelAxesConstraintType2(int*);
    int getTopNumber();
};

#endif
