#ifndef _RBE2_H_
#define _RBE2_H_

#include <Element.d/SuperElement.h>

class RBE2 : public SuperElement
{
  int numcdofs;
  char cdofs[9];

  public:
    RBE2(int numnodes, int *nodenums);

    Element* clone();
    int getTopNumber() { return 503; }
    int numTopNodes() { return nnodes/2; }
    void getGravityForce(CoordSet&, double *, Vector &force, int, GeomState *) { }
    bool isRigidElement() { return true; }

};

#endif

