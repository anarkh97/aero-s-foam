#ifndef _PARALLELAXESCONSTRAINTTYPE1_H_
#define _PARALLELAXESCONSTRAINTTYPE1_H_

#include <Element.d/SuperElement.h>

class ParallelAxesConstraintType1 : public SuperElement
{
  public:
    ParallelAxesConstraintType1(int*);
    int getTopNumber();
};

#endif
