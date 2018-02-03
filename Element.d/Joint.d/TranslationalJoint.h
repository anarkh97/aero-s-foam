#ifndef _TRANSLATIONALJOINT_H_
#define _TRANSLATIONALJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.3
// constrains three rotational dofs

class TranslationalJoint : public SuperElement
{
  public:
    TranslationalJoint(int*);
    int getTopNumber() override;
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
