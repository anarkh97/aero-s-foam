#ifndef _PARALLELAXISCONSTRAINTTYPE1_H_
#define _PARALLELAXISCONSTRAINTTYPE1_H_

#include <Element.d/SuperElement.h>

class ParallelAxisConstraintType1 : public SuperElement
{
  public:
    ParallelAxisConstraintType1(int*);
    int getTopNumber();
};

#endif
