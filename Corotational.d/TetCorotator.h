#ifndef _TET_COROTATOR_H_
#define _TET_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class TetCorotator : public Corotator {
     int nodeNum[4];
     double em;                 // elastic modulus
     double nu;                 // initial cross-sectional area
   public:

     // Constructor
     TetCorotator(int nn[4], double, double, CoordSet &);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f);

     void   formGeometricStiffness(GeomState &gs, CoordSet &cs, 
                                   FullSquareMatrix &elk, double *f);

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     void getNLVonMises(Vector&, Vector&,
                        GeomState &, CoordSet &, int);

     void getNLAllStress(FullM&, Vector&,
                         GeomState &, CoordSet &, int);

     double computeShapeGrad(CoordSet &nodes, double nGrad[4][3]);
 
     void  computeStrainGrad(GeomState &, double nGrad[4][3], double dedU[12][6]);

     void computePiolaStress(GeomState &, CoordSet &cs,
                              double  stress[4][7], double strain[4][7]);

     double getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
