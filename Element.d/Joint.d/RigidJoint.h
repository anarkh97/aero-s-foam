#ifndef _RIGIDJOINT_H_
#define _RIGIDJOINT_H_

#include <Element.d/SuperElement.h>

class RigidJoint : public SuperElement
{
  public:
    RigidJoint(int*);
    int getTopNumber();
};

#endif
