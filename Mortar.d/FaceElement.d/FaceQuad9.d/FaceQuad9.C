// ----------------------------------------------------------------
// HB - 03/01/04
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//          y
//          ^
//          |           Shape functions
//          7           ---------------
//   4+-----+-----+3      Phi[1] = (1/4).((1-x).(1-y).(-1-x-y)+Phi[9])
//    |           |       Phi[2] = (1/4).((1+x).(1-y).(-1+x-y)+Phi[9])
//    |     9     |       Phi[3] = (1/4).((1+x).(1+y).(-1+x+y)+Phi[9])
//   8+     +     +6 -> x Phi[4] = (1/4).((1-x).(1+y).(-1-x+y)+Phi[9])
//    |           |       Phi[5] = (1/2).((1-x.x).(1-y)-Phi[9]) 
//    |           |       Phi[6] = (1/2).((1+x).(1-y.y)-Phi[9]) 
//   1+-----+-----+2      Phi[7] = (1/2).((1-x.x).(1+y)-Phi[9]) 
//          5             Phi[8] = (1/2).((1-x).(1-y.y)-Phi[9]) 
//                        Phi[9] =       (1-x.x).(1-y.y) 
// ----------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// STL
#include <map>

// FEM headers 
#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad9.d/FaceQuad9.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// Extern routine
extern "C" {
void _FORTRAN(qgauss)(int &, int &, int &, int &,
               double &,  double &, double &);
}

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceQuad9::RefCoords[9][2] = {{-1.0,-1.0},
                                     { 1.0,-1.0},
                                     { 1.0, 1.0},
                                     {-1.0, 1.0},
                                     { 0.0,-1.0},
                                     { 1.0, 0.0},
                                     { 0.0, 1.0},
                                     {-1.0, 0.0},
                                     { 0.0, 0.0}};
double* 
FaceQuad9::ViewRefCoords() { return(FaceQuad9::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                              CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad9::FaceQuad9(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
  Nodes[4] = nodenums[4];
  Nodes[5] = nodenums[5];
  Nodes[6] = nodenums[6];
  Nodes[7] = nodenums[7];
  Nodes[8] = nodenums[8];
}

// copy constructor
/*
FaceQuad9::FaceQuad9(const FaceQuad9 &FQ8)
{
  // copy nodes Id
  FQ8.GetNodes(Nodes);
}
*/

// -----------------------------------------------------------------------------------------------------
//                              COPY & CLONE 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
/*
FaceElement* 
FaceQuad9::clone()
{
  return new FaceQuad9(*this);
}
*/

void 
FaceQuad9::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
  Nodes[4] = OldToNewNodeIds[Nodes[4]];
  Nodes[5] = OldToNewNodeIds[Nodes[5]];
  Nodes[6] = OldToNewNodeIds[Nodes[6]];
  Nodes[7] = OldToNewNodeIds[Nodes[7]];
  Nodes[8] = OldToNewNodeIds[Nodes[8]];
}

// -----------------------------------------------------------------------------------------------------
//                              GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF LOCAL METHODS
// -------------------------------
int
FaceQuad9::nQuad4Nodes() { return 4; }

int
FaceQuad9::GetQuad4Node(int i) { return Nodes[i]; } 

void
FaceQuad9::GetQuad4Nodes(int *p, int* renumTable)
{
  if(renumTable){
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
  }
}

void
FaceQuad9::GetQuad4Nodes(int *p, std::map<int,int>&renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceQuad9::nNodes() { return(9); }

int
FaceQuad9::GetNode(int i) { return(Nodes[i]); }

void
FaceQuad9::GetNodes(int *p, int* renumTable)
{
  if(renumTable){
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
    p[4] = renumTable[Nodes[4]];
    p[5] = renumTable[Nodes[5]];
    p[6] = renumTable[Nodes[6]];
    p[7] = renumTable[Nodes[7]];
    p[8] = renumTable[Nodes[8]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
    p[4] = Nodes[4];
    p[5] = Nodes[5];
    p[6] = Nodes[6];
    p[7] = Nodes[7];
    p[8] = Nodes[8];
  }
}

void
FaceQuad9::GetNodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
  p[4] = renumTable[Nodes[4]];
  p[5] = renumTable[Nodes[5]];
  p[6] = renumTable[Nodes[6]];
  p[7] = renumTable[Nodes[7]];
  p[8] = renumTable[Nodes[8]];
}

int 
FaceQuad9::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<9; i++)
    if(gNode==Nodes[i]){ found = true; break; }
  if(!found)  
    filePrint(stderr," *** WARNING: In FaceQuad9::GetNodeIndex: node (%6d) does not belong to this element\n",gNode);
  return(i); 
}

int
FaceQuad9::GetFaceElemType() { return(FaceElement::QUADFACEQ9); }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad9::GetACMEFaceElemType() { return ContactSearch::QUADFACEQ8; }
#else
 int
 FaceQuad9::GetACMEFaceElemType() { return 2; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int 
FaceQuad9::nVertices() { return(nQuad4Nodes()); }

int 
FaceQuad9::GetVertex(int i) { return(GetQuad4Node(i)); }

void 
FaceQuad9::GetVertices(int* p, int* renumTable) { GetQuad4Nodes(p, renumTable); }

void 
FaceQuad9::GetVertices(int* p, std::map<int,int>&renumTable) { GetQuad4Nodes(p, renumTable); }

// As ACME doesn't support the Quad9 face element for
// FaceFaceInteraction (FFI), we will pass to it the Quad4 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Quad9 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad9::GetACMEFFIFaceElemType() { return(ContactSearch::QUADFACEL4); }
#else
 int
 FaceQuad9::GetACMEFFIFaceElemType() { return 1; }
#endif

// -----------------------------------------------------------------------------------------------------
//                          MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad9::GetShapeFct(double *Shape, double *m)
{
  double x = m[0];
  double y = m[1];
  double onequart = 1./4.;
  double onehalf  = 1./2.;

  double xm = 1.-x;
  double xp = 1.+x;
  double ym = 1.-y;
  double yp = 1.+y;
  double xxm= 1.-x*x;
  double yym= 1.-y*y;

  Shape[8] =        xxm*yym;
  Shape[0] = onequart*(xm*ym*(-1.-x-y) + Shape[8]);
  Shape[1] = onequart*(xp*ym*(-1.+x-y) + Shape[8]);
  Shape[2] = onequart*(xp*yp*(-1.+x+y) + Shape[8]);
  Shape[3] = onequart*(xm*yp*(-1.-x+y) + Shape[8]);
  Shape[4] = onehalf*(xxm*ym - Shape[8]);
  Shape[5] = onehalf*(yym*xp - Shape[8]);
  Shape[6] = onehalf*(xxm*yp - Shape[8]);
  Shape[7] = onehalf*(yym*xm - Shape[8]);
}

void
FaceQuad9::GetdShapeFct(double *dShapex, double *dShapey, double *m)
{
  double x = m[0];
  double y = m[1];
  double onequart = 1./4.;
  double onehalf  = 1./2.;

  double xm = 1.-x;
  double xp = 1.+x;
  double ym = 1.-y;
  double yp = 1.+y;
  double xxm= 1.-x*x;
  double yym= 1.-y*y;

  dShapex[8] = -2.*x*yym;  
  dShapex[0] = onequart*(ym*(2.*x+y) + dShapex[8]);
  dShapex[1] = onequart*(ym*(2.*x-y) + dShapex[8]);
  dShapex[2] = onequart*(yp*(2.*x+y) + dShapex[8]);
  dShapex[3] = onequart*(yp*(2.*x-y) + dShapex[8]);
  dShapex[4] = -x*ym        - onehalf*dShapex[8];  
  dShapex[5] =  onehalf*yym - onehalf*dShapex[8];  
  dShapex[6] = -x*yp        - onehalf*dShapex[8];  
  dShapex[7] = -onehalf*yym - onehalf*dShapex[8];  

  dShapey[8] = -2.*y*xxm;  
  dShapey[0] = onequart*(xm*( x+2.*y) + dShapey[8]);
  dShapey[1] = onequart*(xp*(-x+2.*y) + dShapey[8]);
  dShapey[2] = onequart*(xp*( x+2.*y) + dShapey[8]);
  dShapey[3] = onequart*(xm*(-x+2.*y) + dShapey[8]);
  dShapey[4] = -onehalf*xxm - onehalf*dShapey[8];
  dShapey[5] = -y*xp        - onehalf*dShapex[8];  
  dShapey[6] = onehalf*xxm  - onehalf*dShapey[8];
  dShapey[7] = -y*xm        - onehalf*dShapey[8];  
}

double
FaceQuad9::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  GetShapeFct(Shape, m);
  return(GetJacobian(m, cs));
}

void
FaceQuad9::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape fcts derivatives
  double dShapex[9], dShapey[9];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  double X[9], Y[9], Z[9];
  for(int i=0; i<9; i+=3){
    Node &nd0 = cs.getNode(Nodes[i]);
    X[i  ] = nd0.x; Y[i  ] = nd0.y; Z[i  ] = nd0.z;
    Node &nd1 = cs.getNode(Nodes[i+1]);
    X[i+1] = nd1.x; Y[i+1] = nd1.y; Z[i+1] = nd1.z;
    Node &nd2 = cs.getNode(Nodes[i+2]);
    X[i+2] = nd2.x; Y[i+2] = nd2.y; Z[i+2] = nd2.z;
 }

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<9; i+=3){
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2];
    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2];
  }
}

double
FaceQuad9::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // N = dM/dx x dM/dy
  Normal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  Normal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  Normal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  double NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(NormN!=0.0){
   Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
  }
  return(NormN); // TO BE CHECKED ...
}

void
FaceQuad9::GetIsoParamMappingNormalJacobianProduct(double *JNormal, double *m, CoordSet &cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // J*N = dM/dx x dM/dy
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad9::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  double Shape[9];
  GetShapeFct(Shape,m);

  double X[9], Y[9], Z[9];
  for(int i=0; i<9; i+=3){
    Node &nd0 = cs.getNode(Nodes[i]);
    X[i  ] = nd0.x; Y[i  ] = nd0.y; Z[i  ] = nd0.z;
    Node &nd1 = cs.getNode(Nodes[i+1]);
    X[i+1] = nd1.x; Y[i+1] = nd1.y; Z[i+1] = nd1.z;
    Node &nd2 = cs.getNode(Nodes[i+2]);
    X[i+2] = nd2.x; Y[i+2] = nd2.y; Z[i+2] = nd2.z;
  }
                                                                                                                                        
  M[0] = 0.0; M[1] = 0.0; M[2] = 0.0;
  for(int i=0; i<9; i+=3){
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2];
  }
}

void 
FaceQuad9::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

double
FaceQuad9::GetJacobian(double *m, CoordSet &cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // N = dM/dx x dM/dy
  double N[3];
  N[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  N[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  N[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  return(sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}


// -----------------------------------------------------------------------------------------------------
//                              MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceQuad9::ComputeArea(CoordSet &cs,const int ngp=2)
{
  double Area = 0.0;
  double dA, Shape[4];
  double xi, eta, weight, m[2];
	
  for(int i=0;i<ngp;i++){
    for(int j=0;j<ngp;j++){
     	_FORTRAN(qgauss)(ngp,i,ngp,j,&xi,&eta,&weight);
	m[0] = xi; m[1] = eta;
      //dA = ComputeDiffSurfNormaleAndJacobian(Normal, m, cs);  
      dA = GetShapeFctAndJacobian(Shape, m, cs);
    	Area += weight*dA;
    }
  }
  return Area;
} */

// -----------------------------------------------------------------------------------------------------
//                              MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceQuad9::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[9];
  double xi, eta, weight, m[2];

  FullM Mass(8);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<9;i++)
        for(int j=i;j<9;j++) 
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];    
    }
  }
  // lower part by symmetry 
  for(int i=0;i<9;i++)
    for(int j=0;j<i;j++) 
      Mass[i][j] = Mass[j][i];    

  return(Mass);
}

void 
FaceQuad9::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<9; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[9];
  double xi, eta, weight, m[2];

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<9;i++)
        ShapeIntg[i] += rho*weight*dA*Shape[i];
    }
  }
}


// -----------------------------------------------------------------------------------------------------
//                              PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad9::printNodes()
{
  filePrint(stderr,"   # Quad9 face el., nodes = %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
            Nodes[0],Nodes[1],Nodes[2],Nodes[3],Nodes[4],Nodes[5],Nodes[6],Nodes[7],Nodes[8]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad9::print()
{
  printNodes();
}
