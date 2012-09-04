#ifndef _WELDEDJOINT_H_
#define _WELDEDJOINT_H_

#include <Element.d/SuperElement.h>

class WeldedJoint : public SuperElement
{
  public:
    WeldedJoint(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
