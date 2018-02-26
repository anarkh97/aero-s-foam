#ifndef _RIGIDROTNSPRLINK_H_
#define _RIGIDROTNSPRLINK_H_

#include <Element.d/SuperElement.h>

class RigidRotnSprlink : public SuperElement
{
  public:
    RigidRotnSprlink(int*);
    void setProp(StructProp*, bool = false) override;
    int getTopNumber() override { return 122; }
    bool isSafe() const override { return false; }
    bool isSpring() const override { return true; }
    bool isRigidElement() const override { return true; }
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
