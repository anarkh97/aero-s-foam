#ifndef _SHELL3_COROTATOR_H_
#define _SHELL3_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class Node;
class NodeState;

class Shell3Corotator : public Corotator {
   int n1, n2, n3;
   double origK[18][18];
   int fitAlg;
 public:
   Shell3Corotator(int _n1,int _n2,int _n3,FullSquareMatrix &,int fitAlgShell);

   void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *, double dt, double t);

   void getDExternalForceDu(GeomState &geomState, CoordSet &cs,
                                     FullSquareMatrix &elK, double *locF);

   void getInternalForce(GeomState &, CoordSet &, FullSquareMatrix &, double *, double, double);

   void getExternalForce(GeomState &,CoordSet &, double*);
   
   void formGeometricStiffness(GeomState &, CoordSet &, 
                               FullSquareMatrix &, double *);

   double * getOriginalStiffness() { return (double *)origK; } 

   // checked with Haugen's code
   void extractDefDisp(Node &nd1, Node &nd2, Node &nd3, NodeState &ns1,
                       NodeState &ns2, NodeState &ns3,
                       double xl0[3][3], double xln[3][3],
                       double t0[3][3], double t0n[3][3], double vld[18]);
  
   void getGlobalDisp(GeomState& , CoordSet&, Vector& );

   // checked with Haugen's code - but still incorrect
   void formGeometricStiffness(double xl0[3][3],
                          double xln[3][3], double pmat[18][18], 
                          double gmat[3][18], double f[18],
                          double stiffGeo1[18][18], 
                          double stiffGeo2[18][18], double fe[18]);

   // thanks to Joe Pajot
   void formCorrectGeometricStiffness(double rotvar[3][3][3],
                          double xln[3][3], double pmat[18][18], 
                          double gmat[3][18], double f[18],
                          double stiffGeo1[18][18], 
                          double stiffGeo2[18][18], double fe[18],
			  double t0n[3][3]);
   
   // checked SpinAxialAndMoment with Haugen's code
   void spinAxialAndMoment(double f[], double fnm[][3]);

   // checked SpinAxial with Haugen's code
   void spinAxial(double f[], double fn[][3]);
     
   // checked with Haugen's code
   void formRotationGradientMatrix(double xdij[3][3], 
        double ydij[3][3], double xln[3][3], double gmat[3][18]); 

   // checked gradDefDisp with Haugen's code
   void gradDefDisp(double xl0[][3], double xln[][3], 
                    double pmat[18][18], double gmat[3][18]);

   // checked localCoord with Haugen's code
   void localCoord(double x0[3][3], double xn[3][3],
                   double t0[3][3], double t0n[3][3], double xl0[3][3], 
                   double xln[3][3]);

   void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                            int &nlflag);

   void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                               double *vlr);

   void getNLVonMises(Vector&, Vector& weight,
                      GeomState &, CoordSet &, int);

   void getNLAllStress(FullM&, Vector&,
                       GeomState &, CoordSet &, int);

   double getElementEnergy(GeomState &, CoordSet &);
   
   void reBuildorigK(FullSquareMatrix &);

};

#endif
