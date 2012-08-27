#ifndef _DOTTYPE1CFE_H_
#define _DOTTYPE1CFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/DotType1ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType1ConstraintElement : public ConstraintFunctionElement<DotType1ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis1, axis2;
    double d0;

  public:
    DotType1ConstraintElement(int*, int, int, double=0, int=2); 
    ~DotType1ConstraintElement();
    void buildFrame(CoordSet&);
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet&, Eigen::Array<double,7,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
