#ifndef _DOTTYPE3CFE_H_
#define _DOTTYPE3CFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/DotType3ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType3ConstraintElement : public ConstraintFunctionElement<DotType3ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis;
    double d0;

  public:
    DotType3ConstraintElement(int*, int, int=2);
    ~DotType3ConstraintElement();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void setConstantTerm(double _d0) { d0 = _d0; }

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,10,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
