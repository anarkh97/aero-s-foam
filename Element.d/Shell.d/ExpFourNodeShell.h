#ifndef _EXPFOURNODESHELL_H_
#define _EXPFOURNODESHELL_H_

#include <Element.d/SuperElement.h>

class ExpFourNodeShell : public SuperElement 
{
  public:
    ExpFourNodeShell(int *nodenums);

    void buildFrame(CoordSet& cs);
    bool isShell();
    int  getTopNumber();
    bool hasRot() { return true; }

};

#endif
