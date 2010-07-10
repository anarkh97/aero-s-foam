#ifndef _RIGIDSPRING_H_
#define _RIGIDSPRING_H_

#include <Element.d/SuperElement.h>

class RigidSpring : public SuperElement
{
  public:
    RigidSpring(int*);
    int getTopNumber();
    bool isRigidElement() { return true; }
    bool isSpring() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif