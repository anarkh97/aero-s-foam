#ifndef _RIGIDROTNSPRLINK_H_
#define _RIGIDROTNSPRLINK_H_

#include <Element.d/SuperElement.h>

class RigidRotnSprlink : public SuperElement
{
  public:
    RigidRotnSprlink(int*);
    void setProp(StructProp*, bool = false);
    int getTopNumber() override { return 122; }
    bool isSafe() { return false; }
    bool isSpring() { return true; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
