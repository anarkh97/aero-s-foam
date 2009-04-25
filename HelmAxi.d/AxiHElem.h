#ifndef  _AXIHELEM_H_
#define  _AXIHELEM_H_

#include <Element.d/Element.h>
#include <stdio.h>

class AxiHElement: public Element {
  public:
    virtual FullSquareMatrix stiffteta(CoordSet&, double *d) = 0;
    virtual void buildMesh3D(int &eNum, FILE *, int nodeInc, int nSlices) = 0;
    virtual void buildMesh2D(int &eNum, FILE *, int nodeInc, int nSlices) = 0;
};
#endif
