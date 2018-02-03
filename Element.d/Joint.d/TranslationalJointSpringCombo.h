#ifndef _TRANSLATIONALJOINTSPRINGCOMBO_H_
#define _TRANSLATIONALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class TranslationalJointSpringCombo : public SuperElement
{
  public:
    TranslationalJointSpringCombo(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
