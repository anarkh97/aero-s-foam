#ifndef _EIGHTNODEBRICKRIGIDMPC_H_
#define _EIGHTNODEBRICKRIGIDMPC_H_

#include <Element.d/SuperElement.h>

class EightNodeBrickRigidMpc : public SuperElement
{
  public:
    EightNodeBrickRigidMpc(int*);
    int getTopNumber();
};

#endif

