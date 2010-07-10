#ifndef _PRISMATICJOINT_H_
#define _PRISMATICJOINT_H_

#include <Element.d/SuperElement.h>

class PrismaticJoint : public SuperElement
{
  public:
    PrismaticJoint(int*);
    int getTopNumber();
};

#endif