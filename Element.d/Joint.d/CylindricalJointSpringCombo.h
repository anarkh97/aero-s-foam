#ifndef _CYLINDRICALJOINTSPRINGCOMBO_H_
#define _CYLINDRICALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class CylindricalJointSpringCombo : public SuperElement
{
  public:
    CylindricalJointSpringCombo(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
