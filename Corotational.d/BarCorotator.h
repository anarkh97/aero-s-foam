#ifndef _BAR_COROTATOR_H_
#define _BAR_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class BarCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     double em;                 // elastic modulus
     double a0;                 // initial cross-sectional area
     double l0;			// initial length
     double preload;		// preload
   public:

     // Constructor
     BarCorotator(int node1, int node2, double em,
                  double area, double preload, CoordSet &cs);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f);

     void   formInternalForce(double t[3], double p, double *f);

     void   formTangentStiffness(double t[3], double p, double ld, 
                                 double kt[6][6]);

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

     double getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
