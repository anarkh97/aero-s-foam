#ifndef _PRISMATICJOINTSPRINGCOMBO_H_
#define _PRISMATICJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class PrismaticJointSpringCombo : public SuperElement
{
  public:
    PrismaticJointSpringCombo(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
