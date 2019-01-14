#ifndef _PLANARJOINT_H_
#define _PLANARJOINT_H_

#include <Element.d/SuperElement.h>

// constrains one translational and two rotational dofs

class PlanarJoint : public SuperElement
{
public:
	PlanarJoint(int*);
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
