#ifndef _REVOLUTEACTUATOR_H_
#define _REVOLUTEACTUATOR_H_

#include <Element.d/SuperElement.h>

class RevoluteActuator : public SuperElement
{
  public:
    explicit RevoluteActuator(int*);
    void setProp(StructProp *p, bool myProp) override;
    int getTopNumber() override;
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
