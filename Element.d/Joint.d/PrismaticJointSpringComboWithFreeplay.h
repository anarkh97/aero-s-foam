#ifndef _PRISMATICJOINTSPRINGCOMBOWITHFREEPLAY_H_
#define _PRISMATICJOINTSPRINGCOMBOWITHFREEPLAY_H_

#include <Element.d/SuperElement.h>

class PrismaticJointSpringComboWithFreeplay : public SuperElement
{
public:
	PrismaticJointSpringComboWithFreeplay(int*);
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
