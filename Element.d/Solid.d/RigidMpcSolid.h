#ifndef _RIGIDMPCSOLID_H_
#define _RIGIDMPCSOLID_H_

#include <Element.d/SuperElement.h>

class RigidMpcSolid : public SuperElement
{
  public:
    RigidMpcSolid(int numnodes, int *nodenums);

    Element* clone();
    int getTopNumber() { return 101; }
    int numTopNodes() { return 2; }

};

#endif

