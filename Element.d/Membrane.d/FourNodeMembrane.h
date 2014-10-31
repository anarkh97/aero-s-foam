#ifndef _FOURNODEMEMBRANE_H_
#define _FOURNODEMEMBRANE_H_

#include <Element.d/SuperElement.h>

class FourNodeMembrane : public SuperElement 
{
  public:
    FourNodeMembrane(int *nodenums);

    int getTopNumber();
    PrioInfo examine(int sub, MultiFront *mf);
};

#endif
