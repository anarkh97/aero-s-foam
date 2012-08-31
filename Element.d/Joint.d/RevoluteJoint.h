#ifndef _REVOLUTEJOINT_H_
#define _REVOLUTEJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.5
// constrains three translational and two rotational dofs

class RevoluteJoint : public SuperElement
{
  public:
    RevoluteJoint(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
