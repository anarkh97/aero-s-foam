#ifndef _CYLINDRICALJOINT_H_
#define _CYLINDRICALJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.6
// constrains two translational and two rotational dofs

class CylindricalJoint : public SuperElement
{
  public:
    CylindricalJoint(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
