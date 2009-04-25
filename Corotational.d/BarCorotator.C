#include <math.h>

#include <Math.d/FullSquareMatrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/BarCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

BarCorotator::BarCorotator(int _n1, int _n2, double _em, double _a0, 
                           double _preload, CoordSet& cs)
{
 n1 = _n1;	// Node 1
 n2 = _n2;	// Node 2
 em = _em;	// Elastic modulus
 a0 = _a0;	// Original Cross-sectional area
 preload = _preload;  // Preload in Truss

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Compute original length of bar element
 double dx = node2.x - node1.x;
 double dy = node2.y - node1.y;
 double dz = node2.z - node1.z;
 l0 = sqrt(dx*dx + dy*dy + dz*dz);
}

void
BarCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs, 
                               FullSquareMatrix &elK, double *f)
/*******************************************************************
 *
 * Purpose :
 *  Compute Internal force vector and Tangent Stiffness matrix 
 *  for Co-rotated bar element in current configuration.
 *
 * Input Variables:
 *  geomState : current configuration 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * x0[i][j] : global coordinate component j of node i in reference configuration
 * xn[i][j] : global coordinate component j of node i in current configuration
 * t[i]     : transformation matrix (1st vector only) for current state
 * t0[i]    : transformation matrix (1st vector only) for initial state
 * l0       : original length 
 * ld       : deformed length
 * a0       : Cross-sectional area
 * em       : Elastic modulus
 * e        : Green-Lagrange (GL) strain
 * sigma    : PK2 axial stress
 * p        : axial force in local coordinates
 * preload  : axial preload
 * f0       : internal force vector in initial configuration due to preload
 *
 * Output :
 *  f        : internal force vector in current configuration
 *  elK      : element tangent stiffness matrix
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j;
 double kt[6][6], xn[2][3], x0[2][3], t[3], t0[3], f0[6];

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of C0 configuration
 x0[0][0] = node1.x;
 x0[0][1] = node1.y;
 x0[0][2] = node1.z;
 x0[1][0] = node2.x;
 x0[1][1] = node2.y;
 x0[1][2] = node2.z;

 // Set coordinates of Cn configuration 
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Form transformation tensor: t (1st vector only)
 // for current state 
 t[0] = dx/ld;
 t[1] = dy/ld;
 t[2] = dz/ld;

 // Form transformation tensor: t0 (1st vector only)
 // for initial state 
 dx = x0[1][0] - x0[0][0];
 dy = x0[1][1] - x0[0][1];
 dz = x0[1][2] - x0[0][2];
 t0[0] = dx/l0;
 t0[1] = dy/l0;
 t0[2] = dz/l0;

 // Compute current GL-strain
 double e = (ld - l0)/l0;

 // Compute current PK2-stress
 double sigma = em*e;

 // Compute current axial force: p
 double p = sigma*a0;

 // Add Preload
 p += preload;

 // Form current internal force: f
 formInternalForce(t, p, f);

 // Form initial internal force: f0
 formInternalForce(t0, preload, f0);

 for(i=0; i<6; ++i)
   f[i] -= f0[i];

 // Form tangent stiffness matrix
 formTangentStiffness(t, p, ld, kt);

 // Copy tangent stiffness matrix to element K matrix
 for(i=0; i<6; ++i)
   for(j=0; j<6; ++j)
     elK[i][j] = kt[i][j];
}

void
BarCorotator::formInternalForce(double t[3], double p, double *f)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for Co-rotated bar
 *  element in current configuration.
 *
 * Input :
 *  t[i]     : transformation matrix (1st vector only)
 *  p        : axial force in local coordinates
 *
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
  // Compute internal force in local system and store in f
     f[0] =   p;
     f[1] = 0.0;
     f[2] = 0.0;

  // Transform to global coordinate by Fg = T'*Fl and store in f
     // Shortened form, since f[1] = f[2] = 0.0
     f[3] = t[0]*f[0];
     f[4] = t[1]*f[0];
     f[5] = t[2]*f[0];

     f[0] = -f[3];
     f[1] = -f[4];
     f[2] = -f[5];
}

void
BarCorotator::formTangentStiffness(double t[3], double p, 
                                   double ld, double kt[6][6])
/*******************************************************************
 * 
 * Purpose :
 *  Compute tangential stiffness for Co-rotated bar element
 *  in current configuration.
 *
 * Input :
 *  t        : transformation matrix (1st vector only)
 *  p        : axial force in local coordinates
 *  l0       : original length
 *  ld       : deformed length
 *
 * Output :
 *  kt       : tangent stiffness matrix in current configuration
 *
 *****************************************************************/
{
     int i, j;
     double c1;

  // Zero stiffness matrix
     for(i=0; i<6; ++i )
       for(j=0; j<6; ++j )
         kt[i][j] = 0.0;

  // Compute stiffness matrix in local coordinate system 
     c1 = (a0*em)/l0;

  // Transform to global system by: Kt = T'*Kl*T 
     for(i=0; i < 3; ++i)
       for(j=0; j < 3; ++j){ 
        kt[i][j]  = c1*t[i]*t[j];
        kt[i][j] -= (p/ld)*t[i]*t[j]; //HB
       }
     kt[0][0] += p/ld; 
     kt[1][1] += p/ld; 
     kt[2][2] += p/ld; 
     
  // Fill element stiffness matrix
     for(i=0; i<3; i++) {
       for(j=0; j<3; j++) {
          kt[  i][3+j] = -kt[i][j];
          kt[3+i][  j] =  kt[i][3+j];
          kt[3+i][3+j] =  kt[i][j];
       }
     }
}

void
BarCorotator::formGeometricStiffness( GeomState &geomState, CoordSet &,
                               FullSquareMatrix &kg, double *f)
{
/*******************************************************************
 *
 * Purpose :
 *  Compute Geometric stiffness for Co-rotated bar element
 *  in current configuration.
 *
 * Input :
 * geomState: current state of nodes (coordinates, rotation tensors)
 *
 * Local variables:
 *  t        : transformation matrix (1st vector only)
 *  p        : axial force in local coordinates
 *  l0       : undeformed (initial) length of bar
 *  ld       : deformed length
 *  e        : current GL-strain
 *  sigma    : current PK-2 stress in bar
 *  p        : current axial force in bar
 *  a0       : initial area of bar
 *
 * Output :
 *  kg       : geometric tangent stiffness matrix in current configuration
 *
 *****************************************************************/
 // Declare local variables
 int    i, j;
 double xn[2][3];

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Compute current GL-strain
 double e = (ld - l0)/l0;

 // Compute current PK2-stress
 double sigma = em*e;

 // Compute current axial force: p
 double p = sigma*a0;

 // Add Preload
 p += preload;

 // Zero stiffness matrix
 for(i=0; i<6; ++i )
   for(j=0; j<6; ++j )
     kg[i][j] = 0.0;

 kg[0][0] += p/ld;
 kg[1][1] += p/ld;
 kg[2][2] += p/ld;

 // Fill element stiffness matrix
 for(i=0; i<3; i++) {
   for(j=0; j<3; j++) {
     kg[  i][3+j] = -kg[i][j];
     kg[3+i][  j] =  kg[i][3+j];
     kg[3+i][3+j] =  kg[i][j];
   }
 }
}

void
BarCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, 
                                  double *vld, int &nlflag)
{
 // Set Flag to Use Linear Routines for Stress
 nlflag = 1;

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 double xn[2][3], x0[2][3];

 // Set coordinates of C0 configuration
 x0[0][0] = node1.x;
 x0[0][1] = node1.y;
 x0[0][2] = node1.z;
 x0[1][0] = node2.x;
 x0[1][1] = node2.y;
 x0[1][2] = node2.z;

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);
 double delta_l = ld-l0;

 dx = x0[1][0] - x0[0][0];
 dy = x0[1][1] - x0[0][1];
 dz = x0[1][2] - x0[0][2];

 // scale dx, dy, and dz by the initial length
 dx /= l0;
 dy /= l0;
 dz /= l0;

 // Transform to global coordinate by  dXg = T'*dXl
 // Shortened form, since def[1] = def[2] = 0.0
 vld[3] = dx*delta_l/2.0;
 vld[4] = dy*delta_l/2.0;
 vld[5] = dz*delta_l/2.0;
 vld[0] = -vld[3];
 vld[1] = -vld[4];
 vld[2] = -vld[5]; 
}

void
BarCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                     double *vlr)
{
}

void
BarCorotator::getNLVonMises(Vector &stress, Vector &weight,
                            GeomState &geomState, CoordSet &cs,
                            int strInd)
{
 stress.zero();
 weight.zero();
}

void
BarCorotator::getNLAllStress(FullM &stress, Vector &weight,
                             GeomState &geomState, CoordSet &cs,
                             int strInd)
{
 stress.zero();
 weight.zero();
}
 
double
BarCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given State

 // Get original coordinates of bar's nodes
 //Node &node1 = cs.getNode(n1);
 //Node &node2 = cs.getNode(n2);

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 double xn[2][3];
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Compute current GL-strain
 double e = (ld - l0)/l0;

 // Compute stress
 double sigma = em*e;

 // Add Preload????
 sigma += preload/a0;

 // Compute Energy as 1/2 Integral[e*sigma*dV]
 double Energy = 0.5*e*sigma*a0*l0;

 return Energy;
}
