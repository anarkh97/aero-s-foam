#ifndef _TRANSLATIONALJOINTSPRINGCOMBO_H_
#define _TRANSLATIONALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class TranslationalJointSpringCombo : public SuperElement
{
  public:
    TranslationalJointSpringCombo(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
    void setProp(StructProp *p, bool myProp);
};

#endif
