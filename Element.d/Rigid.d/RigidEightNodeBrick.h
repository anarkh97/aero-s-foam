#ifndef _RIGIDEIGHTNODEBRICK_H_
#define _RIGIDEIGHTNODEBRICK_H_

#include <Element.d/SuperElement.h>

class RigidEightNodeBrick : public SuperElement
{
  public:
    RigidEightNodeBrick(int*);
    int getTopNumber() { return 117; }
    bool isRigidElement() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif

