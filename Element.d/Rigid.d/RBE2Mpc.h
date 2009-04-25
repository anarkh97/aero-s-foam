#ifndef _RBE2MPC_H_
#define _RBE2MPC_H_

#include <Element.d/SuperElement.h>

class RBE2Mpc : public SuperElement
{
  int numcdofs;
  char cdofs[9];

  public:
    RBE2Mpc(int numnodes, int *nodenums);

    Element* clone();
    int getTopNumber() { return 101; }
    int numTopNodes() { return 2; }
};

#endif

