#ifndef _TET10_COROTATOR_H_
#define _TET10_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class Tet10Corotator : public Corotator {
     int nodeNum[10];
     double em;                 // elastic modulus
     double nu;                 // Poisson's ratio
     double Tref;               // ambient temperature
     double alpha;              // thermal expansion coefficient
   public:

     // Constructor
     Tet10Corotator(int nn[10], double, double, CoordSet &, double, double);
     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &gs, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void   getInternalForce(GeomState &gs, CoordSet &cs,
                             FullSquareMatrix &elk, double *f, double dt, double t);

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                              int &nlflag);

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                 double *vlr);

     void getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                        GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
                        double ylayer = 0, double zlayer = 0, int avgnum = 0, int measure = -1);

     void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
                         GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
                         int measure = -1);

     double computeShapeGrad(GeomState &nodes, double nGrad[10][3]);

     double computeStrainGrad(GeomState &geomState, CoordSet &, double dedU[30][6],
                              double m[3]);

     void computePiolaStress(GeomState &, CoordSet &cs,
                             double stress[10][7], double strain[10][7]);

     double getElementEnergy(GeomState &gs, CoordSet &cs);

};

#endif
