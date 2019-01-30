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
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
