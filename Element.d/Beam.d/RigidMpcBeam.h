#ifndef _RIGIDMPCBEAM_H_
#define _RIGIDMPCBEAM_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcBeam : public RigidMpcElement
{
    EFrame *elemframe;
  public:
    RigidMpcBeam(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void computePressureForce(CoordSet&, Vector&, GeomState*, int);
    void setFrame(EFrame *_elemframe) { elemframe = _elemframe; }
    int buildFrame(CoordSet&);
  private:
    void getLength(CoordSet&, double&);
    void updTransMatrix(CoordSet&, GeomState*, double t[3][3], double &);
};

#endif
