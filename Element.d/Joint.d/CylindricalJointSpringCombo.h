#ifndef _CYLINDRICALJOINTSPRINGCOMBO_H_
#define _CYLINDRICALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class CylindricalJointSpringCombo : public SuperElement
{
  public:
    CylindricalJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
    void setProp(StructProp *p, bool myProp);
};

#endif
