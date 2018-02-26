#ifndef _RIGIDTRANSSPRLINK_H_
#define _RIGIDTRANSSPRLINK_H_

#include <Element.d/SuperElement.h>

class RigidTransSprlink : public SuperElement
{
  public:
    RigidTransSprlink(int*);
    void setProp(StructProp*, bool = false) override;
    int getTopNumber() override { return 121; }
    bool isSafe() const override { return false; }
    bool isSpring() const override { return true; }
    bool isRigidElement() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
