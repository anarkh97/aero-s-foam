#ifndef _REVOLUTEACTUATOR_H_
#define _REVOLUTEACTUATOR_H_

#include <Element.d/SuperElement.h>

class RevoluteActuator : public SuperElement
{
  public:
    RevoluteActuator(int*);
    void setProp(StructProp *p, bool myProp);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
