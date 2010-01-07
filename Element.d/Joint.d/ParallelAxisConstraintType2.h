#ifndef _PARALLELAXISCONSTRAINTTYPE2_H_
#define _PARALLELAXISCONSTRAINTTYPE2_H_

#include <Element.d/SuperElement.h>

class ParallelAxisConstraintType2 : public SuperElement
{
  public:
    ParallelAxisConstraintType2(int*);
    int getTopNumber();
};

#endif
