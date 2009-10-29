#ifndef _RIGIDTHREENODESHELL_H_
#define _RIGIDTHREENODESHELL_H_

#include <Element.d/SuperElement.h>

class RigidThreeNodeShell : public SuperElement
{
  public:
    RigidThreeNodeShell(int *nodenums);

    Element* clone();
    int getTopNumber() { return 108; }

};

#endif

