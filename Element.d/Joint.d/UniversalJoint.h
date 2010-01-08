#ifndef _UNIVERSALJOINT_H_
#define _UNIVERSALJOINT_H_

#include <Element.d/SuperElement.h>

class UniversalJoint : public SuperElement
{
  public:
    UniversalJoint(int*);
    int getTopNumber();
};

#endif
