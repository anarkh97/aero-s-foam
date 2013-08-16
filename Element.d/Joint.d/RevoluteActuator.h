#ifndef _REVOLUTEACTUATOR_H_
#define _REVOLUTEACTUATOR_H_

#include <Element.d/SuperElement.h>

class RevoluteActuator : public SuperElement
{
  public:
    RevoluteActuator(int*);
    //void setFrame(EFrame *);
    void setProp(StructProp *p, bool myProp);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
