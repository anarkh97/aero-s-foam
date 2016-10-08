#ifndef _REVOLUTEJOINTSPRINGCOMBOWITHFREEPLAY_H_
#define _REVOLUTEJOINTSPRINGCOMBOWITHFREEPLAY_H_

#include <Element.d/SuperElement.h>

class RevoluteJointSpringComboWithFreeplay : public SuperElement
{
  public:
    RevoluteJointSpringComboWithFreeplay(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
