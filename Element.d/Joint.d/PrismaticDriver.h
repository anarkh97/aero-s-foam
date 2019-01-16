#ifndef _PRISMATICDRIVER_H_
#define _PRISMATICDRIVER_H_

#include <Element.d/SuperElement.h>

class PrismaticDriver : public SuperElement
{
public:
	PrismaticDriver(int*);
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
