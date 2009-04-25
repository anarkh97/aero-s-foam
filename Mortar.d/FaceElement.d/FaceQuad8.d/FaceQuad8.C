// ----------------------------------------------------------------
// HB - 03/01/04
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//          y
//          ^
//          |           Shape functions
//          7           ---------------
//   4+-----+-----+3      Phi[1] = (1/4).(1-x).(1-y).(-1-x-y)
//    |           |       Phi[2] = (1/4).(1+x).(1-y).(-1+x-y) 
//    |           |       Phi[3] = (1/4).(1+x).(1+y).(-1+x+y) 
//   8+     x     +6 -> x Phi[4] = (1/4).(1-x).(1+y).(-1-x+y) 
//    |           |       Phi[5] = (1/2).(1-x.x).(1-y) 
//    |           |       Phi[6] = (1/2).(1+x).(1-y.y) 
//   1+-----+-----+2      Phi[7] = (1/2).(1-x.x).(1+y) 
//          5             Phi[8] = (1/2).(1-x).(1-y.y) 
// ----------------------------------------------------------------
// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// STL
#include <map>

// FEM headers 
#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>

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
double FaceQuad8::RefCoords[8][2] = {{-1.0,-1.0},
                                     { 1.0,-1.0},
                                     { 1.0, 1.0},
                                     {-1.0, 1.0},
                                     { 0.0,-1.0},
                                     { 1.0, 0.0},
                                     { 0.0, 1.0},
                                     {-1.0, 0.0}};
double* 
FaceQuad8::ViewRefCoords() { return(FaceQuad8::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                              CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad8::FaceQuad8(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
  Nodes[4] = nodenums[4];
  Nodes[5] = nodenums[5];
  Nodes[6] = nodenums[6];
  Nodes[7] = nodenums[7];
}

// copy constructor
/*
FaceQuad8::FaceQuad8(const FaceQuad8 &FQ8)
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
FaceQuad8::clone()
{
  return new FaceQuad8(*this);
}
*/

void 
FaceQuad8::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
  Nodes[4] = OldToNewNodeIds[Nodes[4]];
  Nodes[5] = OldToNewNodeIds[Nodes[5]];
  Nodes[6] = OldToNewNodeIds[Nodes[6]];
  Nodes[7] = OldToNewNodeIds[Nodes[7]];
}

// -----------------------------------------------------------------------------------------------------
//                              GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF LOCAL METHODS
// -------------------------------
int
FaceQuad8::nQuad4Nodes() { return 4; }

int
FaceQuad8::GetQuad4Node(int i) { return Nodes[i]; } 

void
FaceQuad8::GetQuad4Nodes(int *p, int* renumTable)
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
FaceQuad8::GetQuad4Nodes(int *p, std::map<int,int>&renumTable)
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
FaceQuad8::nNodes() { return 8; }

int
FaceQuad8::GetNode(int i) { return Nodes[i]; }

void
FaceQuad8::GetNodes(int *p, int* renumTable)
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
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
    p[4] = Nodes[4];
    p[5] = Nodes[5];
    p[6] = Nodes[6];
    p[7] = Nodes[7];
  }
}

void
FaceQuad8::GetNodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
  p[4] = renumTable[Nodes[4]];
  p[5] = renumTable[Nodes[5]];
  p[6] = renumTable[Nodes[6]];
  p[7] = renumTable[Nodes[7]];
}

int 
FaceQuad8::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<8; i++)
    if(gNode==Nodes[i]){ found = true; break; }
  if(!found)  
    filePrint(stderr," ### In FaceQuad8::GetNodeIndex(...): node (%6d) does not belong to this element\n",gNode);
  return i; 
}

int
FaceQuad8::GetFaceElemType() { return FaceElement::QUADFACEQ8; }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad8::GetACMEFaceElemType() { return ContactSearch::QUADFACEQ8; }
#else
 int
 FaceQuad8::GetACMEFaceElemType() { return 2; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int 
FaceQuad8::nVertices() { return nQuad4Nodes(); }

int 
FaceQuad8::GetVertex(int i) { return GetQuad4Node(i); }

void 
FaceQuad8::GetVertices(int* p, int* renumTable) { GetQuad4Nodes(p, renumTable); }

void 
FaceQuad8::GetVertices(int* p, std::map<int,int>&renumTable) { GetQuad4Nodes(p, renumTable); }

// As ACME doesn't support the Quad8 face element for
// FaceFaceInteraction (FFI), we will pass to it the Quad4 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Quad8 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad8::GetACMEFFIFaceElemType() { return  ContactSearch::QUADFACEL4; }
#else
 int
 FaceQuad8::GetACMEFFIFaceElemType() { return 1; }
#endif

// -----------------------------------------------------------------------------------------------------
//                          MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad8::GetShapeFct(double *Shape, double *m)
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

  Shape[0] = onequart*xm*ym*(-1.-x-y);
  Shape[1] = onequart*xp*ym*(-1.+x-y);
  Shape[2] = onequart*xp*yp*(-1.+x+y);
  Shape[3] = onequart*xm*yp*(-1.-x+y);
  Shape[4] = onehalf*xxm*ym;
  Shape[5] = onehalf*yym*xp;
  Shape[6] = onehalf*xxm*yp;
  Shape[7] = onehalf*yym*xm;
}

void
FaceQuad8::GetdShapeFct(double *dShapex, double *dShapey, double *m)
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

  dShapex[0] = onequart*ym*(2.*x+y);
  dShapex[1] = onequart*ym*(2.*x-y);
  dShapex[2] = onequart*yp*(2.*x+y);
  dShapex[3] = onequart*yp*(2.*x-y);
  dShapex[4] = -x*ym;  
  dShapex[5] = onehalf*yym;  
  dShapex[6] = -x*yp;  
  dShapex[7] = -onehalf*yym;  

  dShapey[0] = onequart*xm*( x+2.*y);
  dShapey[1] = onequart*xp*(-x+2.*y);
  dShapey[2] = onequart*xp*( x+2.*y);
  dShapey[3] = onequart*xm*(-x+2.*y);
  dShapey[4] = -onehalf*xxm;
  dShapey[5] = -y*xp;  
  dShapey[6] = onehalf*xxm;
  dShapey[7] = -y*xm;  
}

double
FaceQuad8::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  GetShapeFct(Shape, m);
  return(GetJacobian(m, cs));
}

void
FaceQuad8::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape fcts derivatives
  double dShapex[8], dShapey[8];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);
  Node &nd5 = cs.getNode(Nodes[4]);
  Node &nd6 = cs.getNode(Nodes[5]);
  Node &nd7 = cs.getNode(Nodes[6]);
  Node &nd8 = cs.getNode(Nodes[7]);

  double X[8], Y[8], Z[8];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;
  X[4] = nd5.x; Y[4] = nd5.y; Z[4] = nd5.z;
  X[5] = nd6.x; Y[5] = nd6.y; Z[5] = nd6.z;
  X[6] = nd7.x; Y[6] = nd7.y; Z[6] = nd7.z;
  X[7] = nd8.x; Y[7] = nd8.y; Z[7] = nd8.z;

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<8; i+=4){
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2] + dShapex[i+3]*X[i+3];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2] + dShapex[i+3]*Y[i+3];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2] + dShapex[i+3]*Z[i+3];
    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2] + dShapey[i+3]*X[i+3];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2] + dShapey[i+3]*Y[i+3];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2] + dShapey[i+3]*Z[i+3];
  }
}

double
FaceQuad8::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
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
  return(NormN); // !!! TO CHECK !!!
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad8::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);
  Node &nd5 = cs.getNode(Nodes[4]);
  Node &nd6 = cs.getNode(Nodes[5]);
  Node &nd7 = cs.getNode(Nodes[6]);
  Node &nd8 = cs.getNode(Nodes[7]);

  double Shape[8];
  GetShapeFct(Shape,m);

  double X[8], Y[8], Z[8];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;
  X[4] = nd5.x; Y[4] = nd5.y; Z[4] = nd5.z;
  X[5] = nd6.x; Y[5] = nd6.y; Z[5] = nd6.z;
  X[6] = nd7.x; Y[6] = nd7.y; Z[6] = nd7.z;
  X[7] = nd8.x; Y[7] = nd8.y; Z[7] = nd8.z;

  M[0] = 0.0; M[1] = 0.0; M[2] = 0.0;
  for(int i=0; i<8; i+=4){
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2] + Shape[i+3]*X[i+3];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2] + Shape[i+3]*Y[i+3];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2] + Shape[i+3]*Z[i+3];
  }
}

void 
FaceQuad8::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

double
FaceQuad8::GetJacobian(double *m, CoordSet &cs)
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
FaceQuad8::ComputeArea(CoordSet &cs,const int ngp=2)
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
FaceQuad8::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[8];
  double xi, eta, weight, m[2];

  FullM Mass(8);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<8;i++)
        for(int j=i;j<8;j++) 
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];    
    }
  }
  // lower part by symmetry 
  for(int i=0;i<8;i++)
    for(int j=0;j<i;j++) 
      Mass[i][j] = Mass[j][i];    

  return(Mass);
}

void 
FaceQuad8::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<8; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[8];
  double xi, eta, weight, m[2];

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<8;i++)
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
FaceQuad8::printNodes()
{
  filePrint(stderr,"   # Quad8 face el., nodes = %6d %6d %6d %6d %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],Nodes[3],
                                                                                    Nodes[4],Nodes[5],Nodes[6],Nodes[7]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad8::print()
{
  printNodes();
}
