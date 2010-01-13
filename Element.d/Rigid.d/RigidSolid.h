#ifndef _RIGIDSOLID_H_
#define _RIGIDSOLID_H_

#include <Element.d/SuperElement.h>

class RigidSolid : public SuperElement
{
  public:
    RigidSolid(int, int*);
    void buildFrame(CoordSet& cs);
    int getTopNumber();
    int numTopNodes();
    bool isRigidElement() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif

