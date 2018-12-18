#ifndef _REVOLUTEDRIVER_H_
#define _REVOLUTEDRIVER_H_

#include <Element.d/SuperElement.h>

class RevoluteDriver : public SuperElement
{
  public:
    RevoluteDriver(int*);
    int getTopNumber() const override;
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
