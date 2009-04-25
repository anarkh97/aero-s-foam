#ifndef _TRIANGLETHERMAL_COROTATOR_H_
#define _TRIANGLETHERMAL_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class TriangleThermalCorotator : public Corotator {
     int n1; 			// node number 1
     int n2; 			// node number 2
     int n3;                    // node number 3
     double A;                  // area
     double eps;                // emissivity of the body
     double Tr;                 // Temperature of the enclosure receiving the radiation

   public:

     // Constructor
     TriangleThermalCorotator(int node1, int node2, int node3, double area, double epsilon, double Tr, CoordSet &cs);

     double * getOriginalStiffness() { return (double*) 0; }

     void   getStiffAndForce(GeomState &ts, CoordSet &cs, 
                             FullSquareMatrix &elk, double *f);

     void   formInternalForce(double xn[3], double A, double eps, double Tr, double f[3]);

     void   formTangentStiffness(double xn[3], double A, double eps, double kt[3][3]);

};

#endif
