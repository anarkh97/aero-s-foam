#ifndef _RIGIDTRANSSPRLINK_H_
#define _RIGIDTRANSSPRLINK_H_

#include <Element.d/SuperElement.h>

class RigidTransSprlink : public SuperElement
{
  public:
    RigidTransSprlink(int*);
    void setProp(StructProp*, bool = false);
    int getTopNumber() { return 121; }
    bool isSafe() { return false; }
    bool isSpring() { return true; }
    bool isRigidElement() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
