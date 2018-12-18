#ifndef _UNIVERSALJOINTSPRINGCOMBO_H_
#define _UNIVERSALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class UniversalJointSpringCombo : public SuperElement
{
  public:
    UniversalJointSpringCombo(int*);
    int getTopNumber() const override;
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
