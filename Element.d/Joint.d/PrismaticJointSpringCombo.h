#ifndef _PRISMATICJOINTSPRING_H_
#define _PRISMATICJOINTSPRING_H_

#include <Element.d/SuperElement.h>

class PrismaticJointSpringCombo : public SuperElement
{
  public:
    PrismaticJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
