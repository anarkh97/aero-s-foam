#ifndef _CYLINDRICALJOINT_H_
#define _CYLINDRICALJOINT_H_

#include <Element.d/SuperElement.h>

class CylindricalJoint : public SuperElement
{
  public:
    CylindricalJoint(int*);
    int getTopNumber();
};

#endif
