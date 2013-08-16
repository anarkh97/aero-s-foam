#ifndef _FOLLOWERMOMENTELEMENT_H_
#define _FOLLOWERMOMENTELEMENT_H_

#include <Element.d/Function.d/ExternalForce.d/FollowerMomentForceFunction.h>
#include <Element.d/Force.d/ForceFunctionElement.h>

class FollowerMomentElement : public ForceFunctionElement<Simo::FollowerMomentForceFunction>
{
    Eigen::Matrix3d *C0; // initial frame (axes stored row-wise)

  public:
    static const DofSet NODALINPUTDOFS[1];
    static const DofSet NODALOUTPUTDOFS[1];
    FollowerMomentElement(int* _nn); 
    ~FollowerMomentElement();

    void setFrame(EFrame *);
    bool hasRot() { return true; }

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,12,1>& sconst, Eigen::Array<int,0,1>&,
                      GeomState* = NULL, double = 0);
};

#endif