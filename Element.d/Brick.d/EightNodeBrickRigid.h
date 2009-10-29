#ifndef _EIGHTNODEBRICKRIGID_H_
#define _EIGHTNODEBRICKRIGID_H_

#include <Element.d/SuperElement.h>

class EightNodeBrickRigid : public SuperElement
{
  public:
    EightNodeBrickRigid(int *);

    Element* clone();
    int getTopNumber() { return 117; }

};

#endif

