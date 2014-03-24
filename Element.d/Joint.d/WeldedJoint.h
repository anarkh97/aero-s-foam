#ifndef _WELDEDJOINT_H_
#define _WELDEDJOINT_H_

#include <Element.d/SuperElement.h>

class WeldedJoint : public SuperElement
{
    EFrame *elemframe;
    bool myframe;
  public:
    WeldedJoint(int*);
    ~WeldedJoint();
    void setFrame(EFrame *_elemframe);
    void buildFrame(CoordSet&);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
