#ifndef _BARTHERMAL_COROTATOR_H_
#define _BARTHERMAL_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class BarThermalCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     double P;                  // perimeter
     double eps;                // emissivity of the body
     double Tr;                 // temperature of the enclosure receiving the radiation
     double l0;			// initial length

   public:

     BarThermalCorotator(int node1, int node2, double perimeter, double epsilon, double Tr, CoordSet &cs);

     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &ts, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f);

     void   formInternalForce(double t[2], double l0, double P, double eps, double Tr, double *f);

     void   formTangentStiffness(double xn[2],double l0, 
                                 double P, double eps, double kt[2][2]);

};

#endif
