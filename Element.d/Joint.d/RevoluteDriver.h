#ifndef _REVOLUTEDRIVER_H_
#define _REVOLUTEDRIVER_H_

#include <Element.d/SuperElement.h>

class RevoluteDriver : public SuperElement
{
  public:
    RevoluteDriver(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
