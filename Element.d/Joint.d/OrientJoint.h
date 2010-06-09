#ifndef _ORIENTJOINT_H_
#define _ORIENTJOINT_H_

#include <Element.d/SuperElement.h>

// The oorient joint is a two-node element. In this joint, the relative rotational degrees 
// of freedom are fixed while the displacement degrees of freedom are left free.

class OrientJoint : public SuperElement
{
  public:
    OrientJoint(int*);
    int getTopNumber();
};

#endif
