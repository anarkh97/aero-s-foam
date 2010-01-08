#ifndef _PARALLELAXESCONSTRAINTTYPE2_H_
#define _PARALLELAXESCONSTRAINTTYPE2_H_

#include <Element.d/SuperElement.h>

class ParallelAxesConstraintType2 : public SuperElement
{
  public:
    ParallelAxesConstraintType2(int*);
    int getTopNumber();
};

#endif
