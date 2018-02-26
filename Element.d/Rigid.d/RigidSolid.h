#ifndef _RIGIDSOLID_H_
#define _RIGIDSOLID_H_

#include <Element.d/SuperElement.h>

class RigidSolid : public SuperElement
{
  public:
    RigidSolid(int, int*);
    void buildFrame(CoordSet& cs);
    int getTopNumber() override;
    int numTopNodes() override;
    bool isRigidElement() const override { return true; }
    bool isSafe() const override;
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif

