#ifndef _FOLLOWERFORCEELEMENT_H_
#define _FOLLOWERFORCEELEMENT_H_

#include <Element.d/Function.d/ExternalForce.d/FollowerForceFunction.h>
#include <Element.d/Force.d/ForceFunctionElement.h>

class FollowerForceElement : public ForceFunctionElement<Simo::FollowerForceFunction>
{
    Eigen::Matrix3d *C0; // initial frame (axes stored row-wise)

  public:
    static const DofSet NODALINPUTDOFS[1];
    static const DofSet NODALOUTPUTDOFS[1];
    FollowerForceElement(int* _nn); 
    ~FollowerForceElement();

    void setFrame(EFrame *);
    bool hasRot() { return true; }

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,12,1>& sconst, Eigen::Array<int,0,1>&,
                      GeomState* = NULL, double = 0);
};

#endif
