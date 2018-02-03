#ifndef _PRISMATICJOINTSPRINGCOMBOWITHFREEPLAY_H_
#define _PRISMATICJOINTSPRINGCOMBOWITHFREEPLAY_H_

#include <Element.d/SuperElement.h>

class PrismaticJointSpringComboWithFreeplay : public SuperElement
{
  public:
    PrismaticJointSpringComboWithFreeplay(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
