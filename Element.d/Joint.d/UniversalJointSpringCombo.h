#ifndef _UNIVERSALJOINTSPRINGCOMBO_H_
#define _UNIVERSALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class UniversalJointSpringCombo : public SuperElement
{
  public:
    UniversalJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
    void setProp(StructProp *p, bool myProp);
};

#endif
