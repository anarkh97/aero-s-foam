#ifndef _PARALLELAXESCONSTRAINTTYPE1_H_
#define _PARALLELAXESCONSTRAINTTYPE1_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.2
// Constraint Building Block type BB2, also known as a parallel axes constraint
// two constrained rotational dofs

class ParallelAxesConstraintType1 : public SuperElement
{
  public:
    ParallelAxesConstraintType1(int*);
    int getTopNumber();
};

#endif
