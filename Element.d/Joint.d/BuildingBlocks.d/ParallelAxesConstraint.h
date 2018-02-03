#ifndef _PARALLELAXESCONSTRAINT_H_
#define _PARALLELAXESCONSTRAINT_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.2
// Constraint Building Block type BB2, also known as a parallel axes constraint
// two constrained rotational dofs

class ParallelAxesConstraint : public SuperElement
{
  public:
    ParallelAxesConstraint(int*, int=0);
    int getTopNumber() override;
    bool hasRot() { return true; }
};

#endif
