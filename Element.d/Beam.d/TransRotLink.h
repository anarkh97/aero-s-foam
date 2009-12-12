#ifndef _TRANSROTLINK_H_
#define _TRANSROTLINK_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class TransRotLink : public RigidMpcElement
{
    EFrame *elemframe;
  public:
    TransRotLink(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void computePressureForce(CoordSet&, Vector&, GeomState*, int);
    void setFrame(EFrame *_elemframe) { elemframe = _elemframe; }
    int buildFrame(CoordSet&);
    //void getStiffAndForce(GeomState& gState, CoordSet& cs, FullSquareMatrix& Ktan, double* f);
    void updateLMPCs(GeomState& gState, CoordSet& cs);
    void getJacobian(GeomState& gState, CoordSet&, int, FullSquareMatrix& J);
  private:
    void getLength(CoordSet&, double&);
    void updTransMatrix(CoordSet&, GeomState*, double t[3][3], double &);
};

#endif
