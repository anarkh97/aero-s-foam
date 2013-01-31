#ifndef _BRICK_COROTATOR_H_
#define _BRICK_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class BrickCorotator : public Corotator {
     int nodeNum[8];
     double em;                 // elastic modulus
     double nu;                 // initial cross-sectional area
   public:

     // Constructor
     BrickCorotator(int nn[4], double, double, CoordSet &);
     double * getOriginalStiffness() { return (double*) 0; }

     void     getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void     getInternalForce(GeomState &gs, CoordSet &cs,
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void     extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void     extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     void     getNLVonMises(Vector&, Vector& weight,
                        GeomState &, CoordSet &, int);

     void     getNLAllStress(FullM&, Vector&,
                         GeomState &, CoordSet &, int);

     double   computeShapeGrad(GeomState &nodes, double nGrad[8][3]);

     double   computeStrainGrad(GeomState &geomState, CoordSet &, double dedU[24][6],
                              int, int, int);

     void     computePiolaStress(GeomState &, CoordSet &cs,
                              double  stress[8][7], double strain[8][7]);

     double   getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
