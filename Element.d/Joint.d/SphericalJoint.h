#ifndef _SPHERICALJOINT_H_
#define _SPHERICALJOINT_H_

#include <Element.d/SuperElement.h>

class SphericalJoint : public SuperElement
{
  public:
    SphericalJoint(int*);
    int getTopNumber();
};

#endif
