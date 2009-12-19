#ifndef _RIGIDMPCBEAM_H_
#define _RIGIDMPCBEAM_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcBeam : public RigidMpcElement
{
    EFrame *elemframe;
    double l0;       // initial length
    double c0[3][3]; // initial frame (axes stored row-wise)
  public:
    RigidMpcBeam(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void computePressureForce(CoordSet&, Vector&, GeomState*, int);
    void setFrame(EFrame *_elemframe) { elemframe = _elemframe; }
    int buildFrame(CoordSet&);
    //void getStiffAndForce(GeomState& gState, CoordSet& cs, FullSquareMatrix& Ktan, double* f);
    void updateLMPCs(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, int, FullSquareMatrix& H);
  private:
    void getLength(CoordSet&, double&);
    void updTransMatrix(CoordSet&, GeomState*, double t[3][3], double &);
};

#endif
