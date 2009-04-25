#include <stdio.h>
#include <math.h>

#include <Corotational.d/Shell3Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Element.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/linkfc.h>

// Define FORTRAN routines as external functions

extern "C" {

  void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                const int &, const double &, double *, const int &,
                double *, const int &, const double &, double *, const int &);

  void _FORTRAN(dgemv)(const char &, const int &,const int &,
                const double &, double *, const int &,
                double *, const int &, const double &, double *, const int &);

  void _FORTRAN(trirotation)(double *, double *);
                              
                             
}

Shell3Corotator::Shell3Corotator(int _n1, int _n2, int _n3, 
                                 FullSquareMatrix &_origK, int fitAlgShell)
{
 // Set Node Numbers
 n1 = _n1;
 n2 = _n2;
 n3 = _n3;

 // Copy Element's stiffness matrix to origK (original stiffness matrix)
 int i, j;
 for(i=0; i<18; ++i)
   for(j=0; j<18; ++j)
     origK[i][j] = _origK[i][j];

 // Set fit algorithm
 fitAlg = fitAlgShell;
}

void
Shell3Corotator::getStiffAndForce(GeomState &geomState, CoordSet &cs, 
                                  FullSquareMatrix &elK, double *f)
{

 // Get Nodes original coordinates (C0 configuration)
 Node &node1 = cs.getNode( n1 );
 Node &node2 = cs.getNode( n2 );
 Node &node3 = cs.getNode( n3 );

 // Get Nodes current coordinates (C0n configuration)
 NodeState &ns1 = geomState[ n1 ];
 NodeState &ns2 = geomState[ n2 ];
 NodeState &ns3 = geomState[ n3 ];
  
 double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18], locF[18];

 // C0    = initial configuration
 // C0n   = nth configuration
 // xl0   = C0 local coordinates
 // xln   = C0n local coordinates
 // t0    = transformation matrix to C0
 // t0n   = transformation matrix to C0n
 // vld   = local deformation vector
 // locF  = local unprojected internal force
 // origK = original stiffness matrix
 
 // Extract deformational displacement from C0 to C0n configurations

 extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld );

 int i, j;

 // Form unprojected internal forces and initialize stiffness matrix

 for(i=0; i<18; ++i) {
  locF[i] = 0.0;
  for(j=0; j<18; ++j)
    elK[i][j] = origK[i][j];
 }

 // compute locF (local Force) as origK*vld
 _FORTRAN(dgemv)('N',18,18,1.0,(double *)origK,18,vld,1,0.0,locF,1);

 // Compute gradients of the nodal deformational pseudorotations
 // Correct element stiffness and internal force

 double rotvar[3][3][3];

 int inode;
 for(inode=0; inode<3; ++inode)
   pseudorot_var( vld+inode*6+3, rotvar[inode] );

  leftmult_rotvar( 3, 1, rotvar, elK);
 rightmult_rotvar( 3, 0, rotvar, elK);

 double fe[18];

 for(inode=0; inode<3; ++inode)
   for(i=0; i<3; ++i) {
     fe[6*inode+i]   = locF[6*inode+i];
     fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                       rotvar[inode][1][i]*locF[6*inode+4] +
                       rotvar[inode][2][i]*locF[6*inode+5];
   }

 // Add second variation of pseudorotations contracted with the
 // nodal moment to the diagonal blocks of the stiffness

 for(inode = 0; inode<3; ++inode) {
   pseudorot_2var(vld+inode*6+3, locF+inode*6+3, rotvar[inode]);
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       elK[i+inode*6+3][j+inode*6+3] += rotvar[inode][i][j];
 }


 // Compute nonlinear projector matrix relative to deformed element
 // and correct stiffness and force
 
 double pmat[18][18], gmat[3][18];

 gradDefDisp(xl0, xln, pmat, gmat);

 double scrstiff[18][18];

 // Form: [K] = [P'][K]

 _FORTRAN(dgemm)('N','T',18,18,18,1.0,elK.data(),18,
                   (double*)pmat,18,0.0,(double*)scrstiff,18);

 // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',18,18,1.0,(double *)pmat,18,fe,1,0.0,f,1);

 // Form geometric stiffness from internal force and material stiffness

 double stiffGeo1[18][18], stiffGeo2[18][18];

 //
 // For a zero deformation, stiffGeo1 and stiffGeo2 are zero matrics.
 //
 formGeometricStiffness(xl0,xln,pmat,gmat,f,stiffGeo1,stiffGeo2);

 // Assemble element tangent Stiffness Matrix K = Kg1 + Kg2 + P'KP

 // Sum geometric stiffness contributions into elK first

 for(i=0; i<18; ++i)
   for(j=0; j<18; ++j)
     elK[i][j] = stiffGeo1[i][j] + stiffGeo2[i][j];

 //
 // Now multiply scrstiff by pmat on the right and sum into elK
 //

 _FORTRAN(dgemm)('N','N',18,18,18,1.0,(double*)pmat,18,
                 (double*)scrstiff,18,1.0,elK.data(),18);

 // transform stiffness matrix from local to global coordinates

  _FORTRAN(trirotation)( elK.data(), (double*)t0n );

 // transform internal force vector from local to global coordinates

 tran_force(f, t0n, 3);

}

void
Shell3Corotator::formGeometricStiffness(GeomState &geomState, CoordSet &cs,
                                  FullSquareMatrix &elK, double *f)
{

 // Get Nodes original coordinates (C0 configuration)
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);
 Node &node3 = cs.getNode(n3);

 // Get Nodes current coordinates (C0n configuration)
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 NodeState &ns3 = geomState[n3];

 double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18];

 // C0    = initial configuration
 // C0n   = nth configuration
 // xl0   = C0 local coordinates
 // xln   = C0n local coordinates
 // t0    = transformation matrix to C0
 // t0n   = transformation matrix to C0n
 // vld   = local deformation vector
 // origK = original stiffness matrix

 // Extract deformational displacement from C0 to C0n configurations

 extractDefDisp(node1,node2,node3,ns1,ns2,ns3,xl0,xln,t0,t0n,vld);

 // Compute nonlinear projector matrix relative to deformed element
 // and correct stiffness and force

 double pmat[18][18], gmat[3][18];

 gradDefDisp(xl0, xln, pmat, gmat);

 double locF[18], stiffGeo1[18][18], stiffGeo2[18][18];

 // compute locF (local Force) as origK*vld
 _FORTRAN(dgemv)('N',18,18,1.0,(double *)origK,18,vld,1,0.0,locF,1);

 double rotvar[3][3][3];

 int inode;
 for(inode=0; inode<3; ++inode)
   pseudorot_var(vld+inode*6+3, rotvar[inode]);

 double fe[18];

 int i;
 for(inode=0; inode<3; ++inode)
   for(i=0; i<3; ++i) {
     fe[6*inode+i]   = locF[6*inode+i];
     fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                       rotvar[inode][1][i]*locF[6*inode+4] +
                       rotvar[inode][2][i]*locF[6*inode+5];
   }

 // Form: {f} = [P']{fe}

 _FORTRAN(dgemv)('N',18,18,1.0,(double *)pmat,18,fe,1,0.0,f,1);

 //
 // For a zero deformation, stiffGeo1 and stiffGeo2 are zero matrics.
 //

 formGeometricStiffness(xl0,xln,pmat,gmat,f,stiffGeo1,stiffGeo2);

 int j;
 for(i=0; i<18; ++i)
   for(j=0; j<18; ++j)
     elK[i][j] = stiffGeo1[i][j] + stiffGeo2[i][j];

 // transform stiffness matrix from local to global coordinates

  _FORTRAN(trirotation)(elK.data(),(double*)t0n);

}

void
Shell3Corotator::extractDefDisp(Node &nd1, Node &nd2, Node &nd3,NodeState &ns1,
                           NodeState &ns2, NodeState &ns3, 
                           double xl0[3][3], double xln[3][3], 
                           double t0[3][3], double t0n[3][3], double vld[18])
/****************************************************************
 *
 *  Purpose:
 *     Form the deformational displacement vector for an element
 *     by subtracting the rigid body displacements.
 *     The deformational nodal rotations are supposed to have
 *     vector properties.
 * 
 *     routine corot_defdisp in C++ programming
 *  Input:
 *     element : element data structure
 *
 *  Output:
 *     xl0  :  xl0[i][j] is coordinate component j of node i.
 *             Local coordinate system coordinates of initial element
 *     xln  :  xln[i][j] is coordinate component j of node i.
 *             Local coordinate system coordinates of deformed element.
 *     t0   :  t0[i][j] is vector component j of basevector i
 *             for the initial element. i.e.
 *             transformation matrix to global.
 *     t0n  :  t0n[i][j] is vector component j of basevector i
 *             for the corotated and deformed element. i.e.
 *             transformation matrix to global.
 *     vld  :  deformational displacement vector for the element in
 *             local coordinate system.
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
 int i, j, inode;
 double x0[3][3], xn[3][3], (* rot[3])[3][3];

 // Set original configuration in x0
 x0[0][0] = nd1.x; // x coordinate of node 1
 x0[0][1] = nd1.y; // y coordinate of node 1
 x0[0][2] = nd1.z; // z coordinate of node 1

 x0[1][0] = nd2.x; // x coordinate of node 2
 x0[1][1] = nd2.y; // y coordinate of node 2
 x0[1][2] = nd2.z; // z coordinate of node 2

 x0[2][0] = nd3.x; // x coordinate of node 3
 x0[2][1] = nd3.y; // y coordinate of node 3
 x0[2][2] = nd3.z; // z coordinate of node 3

 // Set current configuration in xn
 xn[0][0] = ns1.x; // x coordinate of node state 1
 xn[0][1] = ns1.y; // y coordinate of node state 1
 xn[0][2] = ns1.z; // z coordinate of node state 1

 xn[1][0] = ns2.x; // x coordinate of node state 2
 xn[1][1] = ns2.y; // y coordinate of node state 2
 xn[1][2] = ns2.z; // z coordinate of node state 2

 xn[2][0] = ns3.x; // x coordinate of node state 3
 xn[2][1] = ns3.y; // y coordinate of node state 3
 xn[2][2] = ns3.z; // z coordinate of node state 3

 rot[0]   = &(ns1.R); // rotation tensor of node state 1
 rot[1]   = &(ns2.R); // rotation tensor of node state 2
 rot[2]   = &(ns3.R); // rotation tensor of node state 3
 
 // Compute local coordinates
 localCoord( x0, xn, t0, t0n, xl0, xln );

 // Form the displacement from C0n to Cn
 
 // Translation for node 1
 vld[0]  = xln[0][0] - xl0[0][0];
 vld[1]  = xln[0][1] - xl0[0][1];
 vld[2]  = xln[0][2] - xl0[0][2];

 // Translation for node 2
 vld[6]  = xln[1][0] - xl0[1][0];
 vld[7]  = xln[1][1] - xl0[1][1];
 vld[8]  = xln[1][2] - xl0[1][2];

 // Translation for node 3
 vld[12] = xln[2][0] - xl0[2][0];
 vld[13] = xln[2][1] - xl0[2][1];
 vld[14] = xln[2][2] - xl0[2][2];

 // Create rotation part of the deformation vector
 double rot0[3][3], dr[3][3];

 for(inode = 0; inode < 3; ++inode) {

    for(i=0; i<3; ++i)
      for(j=0; j<3; ++j)
        rot0[i][j] = (*rot[inode])[i][0]*t0[j][0]
                    +(*rot[inode])[i][1]*t0[j][1]
                    +(*rot[inode])[i][2]*t0[j][2];

    for(i=0; i<3; ++i)
      for(j=0; j<3; ++j)
        dr[i][j] = t0n[i][0]*rot0[0][j]
                  +t0n[i][1]*rot0[1][j]
                  +t0n[i][2]*rot0[2][j];

    mat_to_vec(dr, vld + inode*6 + 3);

 }

}
	

void
Shell3Corotator::formGeometricStiffness(double [3][3], 
                            double [3][3], double pmat[18][18], 
                            double gmat[3][18], double f[18], 
                            double stiffGeo1[18][18], 
                            double stiffGeo2[18][18])
/****************************************************************
 *
 *  Purpose:
 *     Form geometric stiffness in local coordinate system for
 *     a corotational element based on current material stiffness
 *     and internal force.
 *
 *  Input:
 *
 *     xl0    : local coordinates of undeformed shadow element
 *              xl0[i][j] is local coordinate component j of node i
 *              element is assumed best fit in xy-plane. i.e.
 *              z coordinate is zero for 3-node element.
 *     xln    : local coordinates of deformed element
 *     pmat   : nonlinear projector matrix for current configuration
 *     pmatold: nonlinear projector matrix, not corrected for pseudovectors
 *     fint   : internal force
 *
 *  Output:
 *     stiffGeo1  : geometric stiffness, rotation contribution
 *     stiffGeo2  : geometric stiffness, equilibrium contribution
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j, k;
   double fspin[18][3], fproj[18][3];

// Fspin with both axial and moment contributions
   spinAxialAndMoment( f, fspin );

// Geometric stiffness contribution Kgeo1 = -Fspin*Gmat

   for( i=0; i<18; ++i ) {
      for( j=0; j<18; ++j ) {
         stiffGeo1[i][j] = -( fspin[i][0]*gmat[0][j]
                             +fspin[i][1]*gmat[1][j]
                             +fspin[i][2]*gmat[2][j] );
      }
   }

// Geometric stiffness contribution Kgeo2 = -Gmat'*Fspin'*Pmat
// where Fspin does not contain any moment contributions

// Fspin with only axial contributions
   spinAxial( f, fspin );

// Compute Fproj' = Fspin'*Pmat

   for( i=0; i<3; ++i ) {
      for( j=0; j<18; ++j ) {
         fproj[j][i] = 0.0;
         for( k=0; k<18; ++k )
            fproj[j][i] += fspin[k][i]*pmat[k][j];
      }
   }

   for( i=0; i<18; ++i )
      for( j=0; j<18; ++j )
         stiffGeo2[i][j] = -( gmat[0][i]*fproj[j][0]
			     +gmat[1][i]*fproj[j][1]
                             +gmat[2][i]*fproj[j][2] );

}

void 
Shell3Corotator::spinAxialAndMoment ( double f[], double fnm[][3])
/****************************************************************
 *
 *  Purpose:
 *     Form Fnm matrix based on internal force vector f
 *     routine form_fnm in c programming
 *
 *  Input:
 *     f      : internal force vector
 *
 *  Output:
 *     fnm    : matrix of | Spin(n) | for each  node ordered as columns
 *                        | Spin(m) |
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int dof, nod;

// Fspin with both axial and moment contributions
   for( nod=0; nod<3; ++nod ) {

   // Zero diagonal elements
      for( dof=0; dof<3; ++dof ) {
         fnm[nod*6    +dof][dof] = 0.0;
         fnm[nod*6 +3 +dof][dof] = 0.0;
      }

   // Nonzero axial force
      fnm[nod*6 +1][2] = -f[nod*6 +0];
      fnm[nod*6   ][2] =  f[nod*6 +1];
      fnm[nod*6   ][1] = -f[nod*6 +2];
      fnm[nod*6 +2][1] =  f[nod*6 +0];
      fnm[nod*6 +2][0] = -f[nod*6 +1];
      fnm[nod*6 +1][0] =  f[nod*6 +2];

   // Nonzero moment force
      fnm[nod*6 +4][2] = -f[nod*6 +3];
      fnm[nod*6 +3][2] =  f[nod*6 +4];
      fnm[nod*6 +3][1] = -f[nod*6 +5];
      fnm[nod*6 +5][1] =  f[nod*6 +3];
      fnm[nod*6 +5][0] = -f[nod*6 +4];
      fnm[nod*6 +4][0] =  f[nod*6 +5];
   }

}


void
Shell3Corotator::formRotationGradientMatrix(double x[3][3], double y[3][3], 
                                            double  [3][3], double gmat[3][18])
/***********************************************************************
 *
 *   Compute rotation gradient matrix for a 3 node element
 *   routine gmat_3node in c programs
 *
 *   Input:
 *     fitalg : == 1 x axis along side 1-2, both for t0 and t0n
 *              == 2 x axis rotated from side 1-2 in so that
 *                   the sum of angles between xl0 and xln side edges
 *                   are zero.
 *              == 3 x axis rotated so that xl0 is rotated relative
 *                   to xln with the continuum mechanics definition
 *                   of rotations
 *   x       :  delta x-coordinates  ( element assumed in xy plane )
 *   y       :  delta y-coordinates  ( delta coordinates of Cn )
 *
 *   Output:
 *   gmat :  rotation gradients of rot_x rot_y rot_z in local system
 *           with respect to the local coord. nodal degree of freedoms
 *           The nodal dof ordering is : tx ty tz rx ry rz for each node
 *
 *   Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 **********************************************************************/
{

   static int p[5]={0,1,2,0,1};
   int    ii, jj, kk, i, j, k;
   double length[3], svx[3], svy[3], area2;

// Initialize all entries in gmat to zero
   for( i=0; i<3; ++i )
     for( j=0; j<18; ++j ) 
       gmat[i][j] = 0.0;

// Compute 2 times the area of the triangle
   area2 = (x[1][0]*y[2][0] - x[2][0]*y[1][0]);

// Compute side lengths, heights and side unit vectors for each node
   for(i=0; i<3; ++i) {
      j = p[i+1];
      k = p[i+2];
      length[i] = sqrt(x[k][j]*x[k][j]+y[k][j]*y[k][j]);
      svx[i]    = x[k][j]/length[i];
      svy[i]    = y[k][j]/length[i];
   }

// Compute the fitalg independent theta_x and theta_y parts of gmat

   for( i=0; i<3; ++i ) {
      j = p[i+1];
      k = p[i+2];
      ii = i*6;
      gmat[0][ii +2] = x[k][j]/area2;
      gmat[1][ii +2] = y[k][j]/area2;
   }

// Compute the fitalg dependent theta_z part of gmat

   if ( fitAlg == 1 ) {             // ********* fitalg == 1 **********
      i = 2;
      j = p[i+1];
      k = p[i+2];
      jj = j*6;
      kk = k*6;
      gmat[2][jj   ] =  svy[i]/length[i];
      gmat[2][jj +1] = -svx[i]/length[i];
      gmat[2][kk   ] = -gmat[2][jj   ];
      gmat[2][kk +1] = -gmat[2][jj +1];
   }

   else if ( fitAlg == 2 ) {         // ******** fitalg == 2 **********
      for( i=0; i<3; ++i ) {
         j = p[i+1];
         k = p[i+2];
         ii = i*6;
         gmat[2][ii   ] = (-svy[j]/length[j] +svy[k]/length[k])/3.0;
         gmat[2][ii +1] = ( svx[j]/length[j] -svx[k]/length[k])/3.0;
      }
   }

   else if ( fitAlg == 3 ) {         // ******** fitalg == 3 **********
      for( i=0; i<3; i++ ) {
         j = p[i+1];
         k = p[i+2];
         ii = i*6;
         gmat[2][ii   ] = x[j][k]/(2.0*area2); 
         gmat[2][ii +1] = y[j][k]/(2.0*area2);
      }
   }

   else {                            // ********* fitalg == 1 == default ***
      i = 2;
      j = p[i+1];
      k = p[i+2];
      jj = j*6;
      kk = k*6;
      gmat[2][jj   ] =  svy[i]/length[i];
      gmat[2][jj +1] = -svx[i]/length[i];
      gmat[2][kk   ] = -gmat[2][jj   ];
      gmat[2][kk +1] = -gmat[2][jj +1];
   }

/*
   int i,j;

// Initialize all entries in gmat to zero
   for( i=0; i<3; i++ )
     for( j=0; j<18; j++ )
       gmat[i][j] = 0.0;

// Compute the area of the triangle
   double area2 = (x[1][0]*y[2][0] - x[2][0]*y[1][0]);

   gmat[0][2]   = x[2][1]/area2;
   gmat[0][8]   = x[0][2]/area2;
   gmat[0][14]  = x[1][0]/area2;

   gmat[1][2]   = y[2][1]/area2;
   gmat[1][8]   = y[0][2]/area2;
   gmat[1][14]  = y[1][0]/area2;

   double len1 = xln[0][0]*xln[0][0]+xln[0][1]*xln[0][1];
   double len2 = xln[1][0]*xln[1][0]+xln[1][1]*xln[1][1];
   double len3 = xln[2][0]*xln[2][0]+xln[2][1]*xln[2][1];

   double coef1 = (xln[0][1] - xln[0][0])/(3.0*len1);
   double coef2 = (xln[1][1] - xln[1][0])/(3.0*len2);
   double coef3 = (xln[2][1] - xln[2][0])/(3.0*len3);

   gmat[2][0]  = coef1;
   gmat[2][1]  = coef1;
   gmat[2][6]  = coef2;
   gmat[2][7]  = coef2;
   gmat[2][12] = coef3;
   gmat[2][13] = coef3;

*/

}


void
Shell3Corotator::spinAxial( double f[], double fn[][3] )
/****************************************************************
 *
 *  Purpose:
 *     Form the Fn matrix based on the internal force vector force
 *
 *  Input:
 *     fint   : internal force vector
 *
 *  Output:
 *     fnm    : matrix of | Spin(n) | for each  node ordered as columns
 *                        |    0    |
 *
 *  Coded by: Bjorn Haugen;  Adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int dof, nod;

// Fspin with only axial contributions

   for( nod=0; nod<3; nod++ ) {

   // Zero diagonal elements
      for( dof=0; dof<3; dof++ ) {
         fn[nod*6    +dof][dof] = 0.0;
         fn[nod*6 +3 +dof][dof] = 0.0;
      }

   // Nonzero axial force
      fn[nod*6 +1][2] = -f[nod*6 +0];
      fn[nod*6   ][2] =  f[nod*6 +1];
      fn[nod*6   ][1] = -f[nod*6 +2];
      fn[nod*6 +2][1] =  f[nod*6 +0];
      fn[nod*6 +2][0] = -f[nod*6 +1];
      fn[nod*6 +1][0] =  f[nod*6 +2];

   // Zero moment force
      fn[nod*6 +4][2] = 0.0;
      fn[nod*6 +3][2] = 0.0;
      fn[nod*6 +3][1] = 0.0;
      fn[nod*6 +5][1] = 0.0;
      fn[nod*6 +5][0] = 0.0;
      fn[nod*6 +4][0] = 0.0;
   }

} 

void
Shell3Corotator::gradDefDisp(double [][3], double xln[][3], 
                             double pmat[18][18], double gmat[3][18])
/***********************************************************************
 *
 *   Compute the gradients of the deformational displacement
 *   vector with respect to the visibel dofs for an element
 *   Nonlinear version
 *
 *   routine pmat_nonlin.c in c programs
 *
 *   Input:
 *     fitalg :  as defined in gmat_3nod and gmat_quad
 *     xl0    :  local coordinates of shadow element C0n
 *     xln    :  local coordinates of deformed element Cn
 *               xl0 and xln are assumed
 *               relative to element centroid and
 *               best fit in x-y coordinate system )
 *   Output:
 *     pmat   :  gradient matrix of deformational displacement vector
 *               with respect to the local coord. nodal degree of freedoms
 *               The nodal dof ordering is : tx ty tz rx ry rz for each node
 *
 *   Coded by Bjorn Haugen;  Adjusted for C++ by Teymour Manzouri
 **********************************************************************/
{
   int    nod, inod, jnod, dof;
   double xnij[3][3], ynij[3][3];

// Compute nodal delta coordinates
   for( inod=0; inod<3; inod++ ) {
      for( jnod=0; jnod<3; jnod++ ) {
         xnij[inod][jnod] = xln[inod][0] - xln[jnod][0];
         ynij[inod][jnod] = xln[inod][1] - xln[jnod][1];
      }
   }

// Gmat( gradient of local rotations) dependent part of projector matrix
   formRotationGradientMatrix( xnij, ynij, xln, gmat);

   for( nod=0; nod<3; nod++ ) {
      for( dof=0; dof<18; dof++ ) {

      // Translation part
         pmat[nod*6  ][dof] = -xln[nod][2]*gmat[1][dof]
                              +xln[nod][1]*gmat[2][dof];
         pmat[nod*6+1][dof] =  xln[nod][2]*gmat[0][dof]
                              -xln[nod][0]*gmat[2][dof];
         pmat[nod*6+2][dof] = -xln[nod][1]*gmat[0][dof]
                              +xln[nod][0]*gmat[1][dof];
      // Rotation part
         pmat[nod*6+3][dof] = -gmat[0][dof];
         pmat[nod*6+4][dof] = -gmat[1][dof];
         pmat[nod*6+5][dof] = -gmat[2][dof];
      }
   }

   // Constant part of projector matrix
   // NOTE: fac = -1.0/number_of_nodes

   double fac = -1.0/3.0;

   for( inod=0; inod<3; inod++ ) {

   // Diagonal translation part
      pmat[inod*6   ][inod*6   ] += 1.0;
      pmat[inod*6 +1][inod*6 +1] += 1.0;
      pmat[inod*6 +2][inod*6 +2] += 1.0;

   // Rotation part
      pmat[inod*6 +3][inod*6 +3] += 1.0;
      pmat[inod*6 +4][inod*6 +4] += 1.0;
      pmat[inod*6 +5][inod*6 +5] += 1.0;

      for( jnod=0; jnod<3; jnod++ ) {

      // Nondiagonal translation part
         pmat[inod*6   ][jnod*6   ] += fac;
         pmat[inod*6 +1][jnod*6 +1] += fac;
         pmat[inod*6 +2][jnod*6 +2] += fac;
      }
   }

}



void
Shell3Corotator::localCoord(double x0[3][3], double xn[3][3], 
                            double t0[3][3], double t0n[3][3], double xl0[3][3],
                            double xln[3][3])
/*****************************************************************
 *
 *  Purpose:
 *     Form the transformation matrix to local coordinate system
 *     for a 3 node element.
 *
 *  Method:
 *     Local z-axis is computed normal to the element both for t0 and t0n.
 *         Local x-axis is long side 1-2 of the x0 coordinates for t0
 *     transformation matrix and xl0 are the coordinates of the element
 *     in this local coordinate system.
 *         Local x-axis is rotated relative to side 1-2 for the deformed
 *     element given by xn coordinates.
 *     xln are the coordinates of the element in this local coordinate system.
 *         The x-axis is rotated so that xl0 are a best fit of the undeformed
 *     element on top of the deformed element given by xln. fitalg
 *     gives which algorithm to use for this rotation.
 *
 *
 *  Input;
 *     fitalg : == 1 x axis along side 1-2, both for t0 and t0n
 *              == 2 x axis rotated from side 1-2 in so that
 *                   the sum of angles between xl0 and xln side edges
 *                   are zero.
 *              == 3 x axis rotated so that xl0 is rotated relative
 *                   to xln with the continuum mechanics definition
 *                   of rotations
 * (conf = configuration)
 *     x0[i][j] : global coordinate componenet j of node i in undeformed conf.
 *     xn[i][j] : global coordinate componenet j of node i in deformed conf.
 *
 *  Output:
 *     t0       : transformation matrix to local undeformed coordinates
 *     t0n      : transformation matrix to local deformed coordinates
 *     xl0[i][j]: local coordinate componenet j of node i in undeformed conf.
 *                also best fit of undef. element on def element.
 *     xln[i][j]: local coordinate componenet j of node i in deformed conf.
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *******************************************************************/
{
   int    i,  nod;
   double xc0[3], xcn[3], alpha;
   double l0, ln, s0[3], sn[3], rten[3][3];

// Compute coordinates of centroid for undeformed and deformed configuration
   for( i=0; i<3; ++i ) {
      xc0[i] = ( x0[0][i] + x0[1][i] + x0[2][i] )/3.0;
      xcn[i] = ( xn[0][i] + xn[1][i] + xn[2][i] )/3.0;
   }

// Compute t0 transformation matrix with x axis along side 1-2 
   for( i=0; i<3; i++ ) t0[0][i] = x0[1][i] - x0[0][i];
   normalize( t0[0] );

// local y axis
   for( i=0; i<3; i++ ) t0[1][i] = x0[2][i] - x0[0][i];
   crossprod( t0[0], t0[1], t0[2] );
   normalize( t0[2] );

// local z axis
   crossprod( t0[2], t0[0], t0[1] );
   normalize( t0[1] );

// Compute local coordinates of undeformed element
   for( nod=0; nod<3; nod++) {
      for( i=0; i<3; i++ ) {
         xl0[nod][i] = t0[i][0]*(x0[nod][0] - xc0[0])
                      +t0[i][1]*(x0[nod][1] - xc0[1])
                      +t0[i][2]*(x0[nod][2] - xc0[2]);
     }
   }

// Compute t0n transformation matrix with x axis along side 1-2
   for( i=0; i<3; i++ ) t0n[0][i] = xn[1][i] - xn[0][i];
   normalize( t0n[0] );

   for( i=0; i<3; i++ ) t0n[1][i] = xn[2][i] - xn[0][i];
   crossprod( t0n[0], t0n[1], t0n[2] );
   normalize( t0n[2] );

   crossprod( t0n[2], t0n[0], t0n[1] );
   normalize( t0n[1] );

// Compute local coordinates of deformed element
   for( nod=0; nod<3; nod++) {
      for( i=0; i<3; i++ ) {
         xln[nod][i] = t0n[i][0]*(xn[nod][0] - xcn[0])
                      +t0n[i][1]*(xn[nod][1] - xcn[1])
                      +t0n[i][2]*(xn[nod][2] - xcn[2]);
     }
   }

   if ( fitAlg == 1 ) {       /*********** fitalg == 1 ************/
      alpha = 0.0;
   }

   else if ( fitAlg == 2 ) {  /*********** fitalg == 2 ************/
   /* Compute angle between side 2-3 for deformed and undeformed */
      for( i=0; i<2; i++ ) {
         s0[i] = xl0[2][i] - xl0[1][i];
         sn[i] = xln[2][i] - xln[1][i];
      }
      l0 = sqrt ( s0[0]*s0[0] + s0[1]*s0[1] );
      ln = sqrt ( sn[0]*sn[0] + sn[1]*sn[1] );
      for( i=0; i<2; i++ ) {
         s0[i] /= l0;
         sn[i] /= ln;
      }
      alpha = asin( sn[0]*s0[1] - s0[0]*sn[1] );

   /* Compute angle between side 3-1 for deformed and undeformed */
      for( i=0; i<2; i++ ) {
         s0[i] = xl0[0][i] -xl0[2][i];
         sn[i] = xln[0][i] -xln[2][i];
      }
      l0 = sqrt ( s0[0]*s0[0] + s0[1]*s0[1] );
      ln = sqrt ( sn[0]*sn[0] + sn[1]*sn[1] );
      for( i=0; i<2; i++ ) {
         s0[i] /= l0;
         sn[i] /= ln;
      }
      alpha += asin( sn[0]*s0[1] - s0[0]*sn[1] );
      alpha = -alpha/3.0;
   }

   else if ( fitAlg == 3 ) {  /*********** fitalg == 3 ************/
      alpha = fitalg3_3nodnew( x0, xn );
   }

   else 
     alpha = 0.0;

// Compute rotation vector and rotation tensor for t0n matrix
   for( i=0; i<3; i++ )
     sn[i] = alpha*t0n[2][i];

   vec_to_mat( sn, rten );

// Rotate x-axis angle alpha about z-axis for t0n
   for( i=0; i<3; i++ ) sn[i] = rten[i][0]*t0n[0][0]
                               +rten[i][1]*t0n[0][1]
                               +rten[i][2]*t0n[0][2];

   for( i=0; i<3; i++ ) 
     t0n[0][i] = sn[i];
     
   normalize( t0n[0] );

// Compute y-axis as cross product of z-axis and x-axis
   crossprod( t0n[2], t0n[0], t0n[1] );
   normalize( t0n[1] );

// Recompute local coordinates for deformed configuration
   for( nod=0; nod<3; nod++) {
      for( i=0; i<3; i++ ) {
         xln[nod][i] = t0n[i][0]*(xn[nod][0] - xcn[0])
                      +t0n[i][1]*(xn[nod][1] - xcn[1])
                      +t0n[i][2]*(xn[nod][2] - xcn[2]);
     }
   }

}

void
Shell3Corotator::extractDeformations( GeomState &geomState, CoordSet &cs, 
                 double *vld, int &nlflag)
{
 // Set Flag to Use Linear Routines for Stress
 nlflag = 1;

 // Get Nodes original coordinates (C0 configuration)
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);
 Node &node3 = cs.getNode(n3);

 // Get Nodes current coordinates (C0n configuration)
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 NodeState &ns3 = geomState[n3];
/*
 cerr << "n1 = " << n1 << ", n2 = " << n2 << ", n3 = " << n3 << endl;
 cerr << "geomState = " << ns1.x << "," << ns1.y << "," << ns1.z << "  " 
      << ns2.x << "," << ns2.y << "," << ns2.z << "  "
      << ns3.x << "," << ns3.y << "," << ns3.z << endl;
*/ 
 double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3];

 // C0    = initial configuration
 // C0n   = nth configuration
 // xl0   = C0 local coordinates
 // xln   = C0n local coordinates
 // t0    = transformation matrix to C0
 // t0n   = transformation matrix to C0n
 // vld   = local deformation vector
 // locF  = local unprojected internal force
 // origK = original stiffness matrix

 // Extract deformational displacement from C0 to C0n configurations

 extractDefDisp(node1,node2,node3,ns1,ns2,ns3,xl0,xln,t0,t0n,vld);

 // transform element displacement vector from local to global coordinates

 tran_force( vld, t0, 3 );

// cerr << "vld = "; for(int i=0; i<18; ++i) cerr << vld[i] << " "; cerr << endl;
}

void
Shell3Corotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                 double *vlr)
{
 cerr << "WARNING: Shell3Corotator::extractRigidBodyMotion(...) is not implemented\n";
/*
 // Get Nodes original coordinates (C0 configuration)
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);
 Node &node3 = cs.getNode(n3);

 // Get Nodes current coordinates (C0n configuration)
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 NodeState &ns3 = geomState[n3];
*/
}
void
Shell3Corotator::getNLVonMises(Vector& stress,Vector& weight,
                               GeomState &geomState, CoordSet &cs,
                               int strInd)
{
 stress.zero();
 weight.zero();
}

void
Shell3Corotator::getNLAllStress(FullM& stress,Vector& weight,
                                GeomState &geomState, CoordSet &cs,
                                int strInd)
{
 stress.zero();
 weight.zero();
}


double
Shell3Corotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given Stat

 // Get Nodes original coordinates (C0 configuration)
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);
 Node &node3 = cs.getNode(n3);

 // Get Nodes current coordinates (C0n configuration)
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 NodeState &ns3 = geomState[n3];

 double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3];
 double vld[18], tmp[18];

 // C0    = initial configuration
 // C0n   = nth configuration
 // xl0   = C0 local coordinates
 // xln   = C0n local coordinates
 // t0    = transformation matrix to C0
 // t0n   = transformation matrix to C0n
 // vld   = local deformation vector
 // locF  = local unprojected internal force
 // origK = original stiffness matrix

 // Extract deformational displacement from C0 to C0n configurations

 extractDefDisp(node1,node2,node3,ns1,ns2,ns3,xl0,xln,t0,t0n,vld);

 // transform element displacement vector from local to global coordinates
 tran_force( vld, t0, 3 );

 // Multiply by the original stiffness matrix
 int i,j;
 for(i=0; i<18; ++i) {
   tmp[i] = 0.0;
   for(j=0; j<18; ++j) {
     tmp[i] += origK[i][j]*vld[j];
   }
 }

 // Compute Energy
 double Energy = 0.0;
 for(i=0; i<18; ++i)
   Energy += tmp[i]*vld[i];

 Energy *= 0.5;
 return Energy;
}
