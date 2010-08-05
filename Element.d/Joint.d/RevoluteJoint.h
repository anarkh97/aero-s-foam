#ifndef _REVOLUTEJOINT_H_
#define _REVOLUTEJOINT_H_

#include <Element.d/SuperElement.h>

class RevoluteJoint : public SuperElement
{
  public:
    RevoluteJoint(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
