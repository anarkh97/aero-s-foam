// ----------------------------------------------------------------
// HB - 08/15/03
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//    s
//    ^
//    |
//   2+
//    |\             Shape functions
//    | \            ---------------
//    |  \             Phi[1] = r*(2.r-1)
//   5+   +4           Phi[2] = s*(2.s-1)
//    |    \           Phi[3] = t*(2.t-1), t = 1-r-s
//    |     \          Phi[4] = 4.r.s
//    |      \         Phi[5] = 4.s.t
//   3+---+---+1 -> r  Phi[6] = 4.t.r
//        6
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers 
#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri6::RefCoords[6][2] = {{ 1.0 , 0.0 },
                                    { 0.0 , 1.0 },
                                    { 0.0 , 0.0 },
                                    { 0.5 , 0.5 },
                                    { 0.0 , 0.5 },
                                    { 0.5 , 0.0 }};
double* 
FaceTri6::ViewRefCoords() { return(FaceTri6::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri6::FaceTri6(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
  Nodes[4] = nodenums[4];
  Nodes[5] = nodenums[5];
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
  Nodes[4] = OldToNewNodeIds[Nodes[4]];
  Nodes[5] = OldToNewNodeIds[Nodes[5]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
int
FaceTri6::nTri3Nodes() { return 3; }

int
FaceTri6::GetTri3Node(int i) { return Nodes[i]; }

void
FaceTri6::GetTri3Nodes(int *p, int* renumTable)
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
  }
}

void
FaceTri6::GetTri3Nodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri6::nNodes() { return 6; }

int
FaceTri6::GetNode(int i) { return Nodes[i]; }

void
FaceTri6::GetNodes(int *p, int* renumTable)
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
    p[4] = renumTable[Nodes[4]];
    p[5] = renumTable[Nodes[5]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
    p[4] = Nodes[4];
    p[5] = Nodes[5];
  }
}

void
FaceTri6::GetNodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
  p[4] = renumTable[Nodes[4]];
  p[5] = renumTable[Nodes[5]];
}

int 
FaceTri6::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<6; i++)
    if(gNode==Nodes[i]){ found = true; break; }
  if(!found)  
    filePrint(stderr," ### In FaceTri3::GetNodeIndex(...): node (%6d) does not belong to this element\n",gNode);
  return(i); 
}

int
FaceTri6::nVertices() { return nTri3Nodes(); }

int
FaceTri6::GetVertex(int i) { return GetTri3Node(i); }

void
FaceTri6::GetVertices(int *p, int* renumTable) { GetTri3Nodes(p, renumTable); }

void
FaceTri6::GetVertices(int *p, std::map<int,int>& renumTable) { GetTri3Nodes(p, renumTable); }

int
FaceTri6::GetFaceElemType() { return FaceElement::TRIFACEQ6; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri6::GetACMEFaceElemType() { return ContactSearch::TRIFACEQ6; }
#else
int
FaceTri6::GetACMEFaceElemType() { return 4; }
#endif

// As ACME doesn't support the Tri6 face element for
// FaceFaceInteraction (FFI), we will pass to it the Tri3 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Tri6 face element has straight edges, but its is an 
// APPROXIMATION in the general case (i.e. curved edges/face). 
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri6::GetACMEFFIFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int 
FaceTri6::GetACMEFFIFaceElemType() { return 3; }
#endif
// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void 
FaceTri6::GetShapeFct(double *Shape, double *m)
{
  double r = m[0];
  double s = m[1];
  double t = 1.-r-s; 

  //Shape[0] = t*(2.*t-1.);
  //Shape[1] = r*(2.*r-1.);
  //Shape[2] = s*(2.*s-1.);
  //Shape[3] = 4.*r*t;
  //Shape[4] = 4.*r*s;
  //Shape[5] = 4.*s*t;
  // !! idem ACME !!
  Shape[0] = r*(2.*r-1.);
  Shape[1] = s*(2.*s-1.);
  Shape[2] = t*(2.*t-1.);
  Shape[3] = 4.*r*s;
  Shape[4] = 4.*s*t;
  Shape[5] = 4.*t*r;
}

void
FaceTri6::GetdShapeFct(double *dShapex, double *dShapey, double *m)
{
  double r = m[0];
  double s = m[1];
  double t = 1.-r-s;

  //dShapex[0] = -4.*t+1. ;
  //dShapex[1] =  4.*r-1. ;
  //dShapex[2] =  0.      ;
  //dShapex[3] =  4.*(t-r);
  //dShapex[4] =  4.*s    ;
  //dShapex[5] = -4.*s    ;
 
  //dShapey[0] = -4.*t+1. ;
  //dShapey[1] =  0.      ;
  //dShapey[2] =  4.*s-1. ;
  //dShapey[3] = -4.*r    ;
  //dShapey[4] =  4.*r    ;
  //dShapey[5] =  4.*(t-s);
  // !! idem ACME !!
  dShapex[0] =  4.*r-1. ;
  dShapex[1] =  0.      ;
  dShapex[2] = -4.*t+1. ;
  dShapex[3] =  4.*s    ;
  dShapex[4] = -4.*s    ;
  dShapex[5] =  4.*(t-r);

  dShapey[0] =  0.      ;
  dShapey[1] =  4.*s-1. ;
  dShapey[2] = -4.*t+1. ;
  dShapey[3] =  4.*r    ;
  dShapey[4] =  4.*(t-s);
  dShapey[5] = -4.*r    ;
}

double
FaceTri6::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  GetShapeFct(Shape, m);  
  return(GetJacobian(m, cs));
}

void
FaceTri6::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape fcts derivatives
  double dShapex[6], dShapey[6];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);
  Node &nd5 = cs.getNode(Nodes[4]);
  Node &nd6 = cs.getNode(Nodes[5]);

  double X[6], Y[6], Z[6];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;
  X[4] = nd5.x; Y[4] = nd5.y; Z[4] = nd5.z;
  X[5] = nd6.x; Y[5] = nd6.y; Z[5] = nd6.z;

  dMdx[0] = 0.0; dMdx[1] = 0.0; dMdx[2] = 0.0;
  dMdy[0] = 0.0; dMdy[1] = 0.0; dMdy[2] = 0.0;
  for(int i=0; i<6; i+=3){
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2];
    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2];
  }
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);
  Node &nd5 = cs.getNode(Nodes[4]);
  Node &nd6 = cs.getNode(Nodes[5]);

  double Shape[6];
  GetShapeFct(Shape,m);

  double X[6], Y[6], Z[6];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;
  X[4] = nd5.x; Y[4] = nd5.y; Z[4] = nd5.z;
  X[5] = nd6.x; Y[5] = nd6.y; Z[5] = nd6.z;

  M[0] = 0.0; M[1] = 0.0; M[2] = 0.0;
  for(int i=0; i<6; i+=3){
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2];
  }
}

void 
FaceTri6::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

double
FaceTri6::GetJacobian(double *m, CoordSet &cs)
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

double
FaceTri6::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
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
  return(NormN);
}

void
FaceTri6::GetIsoParamMappingNormalJacobianProduct(double *JNormal, double *m, CoordSet &cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // J*N = dM/dx x dM/dy
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceTri6::ComputeArea(CoordSet &cs,const int ngp=0)
{
  return GetJacobian(cs);
} 
*/

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceTri6::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass;
  Mass.zero();
  fprintf(stderr," *** WARNING:  FaceTri6::ScalarMass(...): NOT IMPLEMENTED !!!\n");
  return Mass;
}

void 
FaceTri6::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  fprintf(stderr," *** WARNING: FaceTri6::IntegrateShapeFcts(...): NOT IMPLEMENTED !!!\n");
  for(int i=0;i<6;i++) { ShapeIntg[i] = 0.0; }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceTri6::printNodes()
{
  filePrint(stderr,"   # Tri6 face el., nodes = %6d %6d %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],
                                                                           Nodes[3],Nodes[4],Nodes[5]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::print()
{
  printNodes();
}
