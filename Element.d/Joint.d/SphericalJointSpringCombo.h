#ifndef _SPHERICALJOINTSPRINGCOMBO_H_
#define _SPHERICALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class SphericalJointSpringCombo : public SuperElement
{
  public:
    SphericalJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
    void setProp(StructProp *p, bool myProp);
};

#endif
