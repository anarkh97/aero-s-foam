#ifndef _REVOLUTEACTUATOR_H_
#define _REVOLUTEACTUATOR_H_

#include <Element.d/SuperElement.h>

class RevoluteActuator : public SuperElement
{
  public:
    RevoluteActuator(int*);
    int getTopNumber();
};

#endif
