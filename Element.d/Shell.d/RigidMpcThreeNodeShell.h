#ifndef _RIGIDMPCTHREENODESHELL_H_
#define _RIGIDMPCTHREENODESHELL_H_

#include <Element.d/SuperElement.h>

class RigidMpcThreeNodeShell : public SuperElement
{
  public:
    RigidMpcThreeNodeShell(int*);
    int getTopNumber();
};

#endif

