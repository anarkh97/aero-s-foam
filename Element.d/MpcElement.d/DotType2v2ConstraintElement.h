#ifndef _DOTTYPE2V2CFE_H_
#define _DOTTYPE2V2CFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/DotType2v2ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType2v2ConstraintElement : public ConstraintFunctionElement<DotType2v2ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis;

  public:
    DotType2v2ConstraintElement(int*, int, int=2);
    ~DotType2v2ConstraintElement();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,9,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
