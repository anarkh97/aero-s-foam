#ifndef _PININSLOTJOINTSPRINGCOMBO_H_
#define _PININSLOTJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class PinInSlotJointSpringCombo : public SuperElement
{
  public:
    PinInSlotJointSpringCombo(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
