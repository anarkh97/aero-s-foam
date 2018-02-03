#ifndef _PININSLOTJOINT_H_
#define _PININSLOTJOINT_H_

#include <Element.d/SuperElement.h>

// constrains two translational and two rotational dofs
// but unlike the cylindrical joint the free axes of translation
// and rotation are not the same, in fact they are orthogonal

class PinInSlotJoint : public SuperElement
{
  public:
    PinInSlotJoint(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
