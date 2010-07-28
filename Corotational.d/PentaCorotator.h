#ifndef _PENTA_COROTATOR_H_
#define _PENTA_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class PentaCorotator : public Corotator {
     int nodeNum[6];
     double em;                 // elastic modulus
     double nu;                 // initial cross-sectional area
   public:

     // Constructor
     PentaCorotator(int nn[6], double, double, CoordSet &);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt);

     void   formGeometricStiffness(GeomState &gs, CoordSet &cs, 
                                   FullSquareMatrix &elk, double *f);

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     void getNLVonMises(Vector&, Vector& weight,
                        GeomState &, CoordSet &, int);

     void getNLAllStress(FullM&, Vector&,
                         GeomState &, CoordSet &, int);

     double computeShapeGrad(GeomState &nodes, double nGrad[6][3]);

     double computeStrainGrad(GeomState &geomState, CoordSet &, double dedU[18][6],
                              double m[3]);

     void computePiolaStress(GeomState &, CoordSet &cs,
                              double  stress[6][7], double strain[6][7]);

     double getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
