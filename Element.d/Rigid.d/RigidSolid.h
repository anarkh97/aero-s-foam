#ifndef _RIGIDSOLID_H_
#define _RIGIDSOLID_H_

#include <Element.d/SuperElement.h>

class RigidSolid : public SuperElement
{
public:
	explicit RigidSolid(int, int*);
	void buildFrame(CoordSet& cs) override;
	int getTopNumber() const override;
	int numTopNodes() const override;
	bool isRigidElement() const override { return true; }
	bool isSafe() const override;
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif

