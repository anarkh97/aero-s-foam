#ifndef _REVOLUTEJOINTSPRINGCOMBO_H_
#define _REVOLUTEJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class RevoluteJointSpringCombo : public SuperElement
{
  public:
    RevoluteJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
