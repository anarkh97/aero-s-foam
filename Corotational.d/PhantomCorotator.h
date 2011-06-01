#ifndef _PHANTOMCOROTATOR_H_
#define _PHANTOMCOROTATOR_H_

#include <Corotational.d/Corotator.h>

class PhantomCorotator : public Corotator {
  public:
    void getStiffAndForce(GeomState &, CoordSet &, 
                          FullSquareMatrix &, double *, double, double);

    void formGeometricStiffness(GeomState &, CoordSet &, 
                                FullSquareMatrix &, double *); 

    double* getOriginalStiffness();

    void extractDeformations(GeomState &geomState, CoordSet &cs, 
                             double *vld, int &nlflag);
    void extractDeformations(GeomState &geomState, CoordSet &cs,
                             DComplex *vld, int &nlflag);

    void getNLVonMises(Vector&, Vector& weight,
                       GeomState &, CoordSet &, int);
    void getNLVonMises(ComplexVector&, Vector& weight,
                       GeomState &, CoordSet &, int);

    void getNLAllStress(FullM&, Vector&,
                        GeomState &, CoordSet &, int);

    double getElementEnergy(GeomState &, CoordSet &);
 
    void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                double *vlr);
};

#endif