#ifndef _PRISMATICACTUATOR_H_
#define _PRISMATICACTUATOR_H_

#include <Element.d/SuperElement.h>

class PrismaticActuator : public SuperElement
{
  public:
    PrismaticActuator(int*);
    void setProp(StructProp *p, bool myProp);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
