#ifndef _TWONODETRUSSRIGIDMPC_H_
#define _TWONODETRUSSRIGIDMPC_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class TwoNodeTrussRigidMpc : public RigidMpcElement
{
    EFrame *elemframe;
  public:
    TwoNodeTrussRigidMpc(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    bool isSafe();
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
