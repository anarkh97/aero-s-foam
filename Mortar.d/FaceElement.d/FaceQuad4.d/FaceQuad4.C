// ----------------------------------------------------------------
// HB - 05/06/03
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
#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// Extern routine
extern "C" {
void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);
}

#define NICER_IMPLEMENTATION

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceQuad4::RefCoords[4][2] = {{-1.0,-1.0},
                                     { 1.0,-1.0},
                                     { 1.0, 1.0},
                                     {-1.0, 1.0}};

double*
FaceQuad4::ViewRefCoords() { return(FaceQuad4::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad4::FaceQuad4(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
}

// copy constructor
/*
FaceQuad4::FaceQuad4(const FaceQuad4 &FQ4)
{
  // copy nodes Id
  FQ4.GetNodes(Nodes);
}
*/

// -----------------------------------------------------------------------------------------------------
//                                            COPY & CLONE 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
/*
FaceElement* 
FaceQuad4::clone()
{
  return new FaceQuad4(*this);
}
*/

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE  METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void 
FaceQuad4::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceQuad4::nNodes() { return 4; }

int
FaceQuad4::GetNode(int i) { return Nodes[i]; }

void
FaceQuad4::GetNodes(int *p, int* renumTable)
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
FaceQuad4::GetNodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
}

int 
FaceQuad4::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<4; i++)
    if(gNode==Nodes[i]){ found = true; break; }
   if(!found){  
    filePrint(stderr," *** WARNING: FaceQuad4::GetNodeIndex(): node (%6d) does not belong to this element\n",gNode);
    printNodes();
   }
  return i; 
}

int
FaceQuad4::GetFaceElemType() { return FaceElement::QUADFACEL4; }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad4::GetACMEFaceElemType() { return ContactSearch::QUADFACEL4; }
#else
 int
 FaceQuad4::GetACMEFaceElemType() { return 1; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int 
FaceQuad4::nVertices() { return nNodes(); }

int 
FaceQuad4::GetVertex(int i) { return GetNode(i); }

void 
FaceQuad4::GetVertices(int* p, int* renumTable) { GetNodes(p, renumTable); }

void 
FaceQuad4::GetVertices(int* p, std::map<int,int>& renumTable) { GetNodes(p, renumTable); }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad4::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#else
 int
 FaceQuad4::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad4::GetShapeFct(double *Shape, double *m)
{
  double x = m[0];
  double y = m[1];

  double d1 = 0.5*(1.0+x);
  double d2 = 0.5*(1.0+y);
  double d3 = 1.0-d1;
  double d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;
}

void
FaceQuad4::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  double x = m[0];
  double y = m[1];
  double onequart = 1./4.;

  double xm = 1.-x;
  double xp = 1.+x;
  double ym = 1.-y;
  double yp = 1.+y;

  dShapex[0] = -onequart*ym;
  dShapex[1] =  onequart*ym;
  dShapex[2] =  onequart*yp;
  dShapex[3] = -onequart*yp;

  dShapey[0] = -onequart*xm;
  dShapey[1] = -onequart*xp;
  dShapey[2] =  onequart*xp;
  dShapey[3] =  onequart*xm;
}

double
FaceQuad4::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  double x = m[0];
  double y = m[1];

  double d1 = 0.5*(1.0+x);
  double d2 = 0.5*(1.0+y);
  double d3 = 1.0-d1;
  double d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;

  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  double a[4], b[4], c[4];
  a[0] = (Y[1]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[1]-Z[0]);
  a[1] = (Y[1]-Y[0])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[1]-Z[0]);
  a[2] = (Y[2]-Y[3])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[2]-Z[3]);
  a[3] = (Y[2]-Y[3])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[3]);

  b[0] = (Z[1]-Z[0])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[1]-X[0]);
  b[1] = (Z[1]-Z[0])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[1]-X[0]);
  b[2] = (Z[2]-Z[3])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[2]-X[3]);
  b[3] = (Z[2]-Z[3])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[2]-X[3]);

  c[0] = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]);
  c[1] = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]);
  c[2] = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]);
  c[3] = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]);

  double N[3];
  N[0] = Shape[0]*a[0]+Shape[1]*a[1]+Shape[2]*a[2]+Shape[3]*a[3];
  N[1] = Shape[0]*b[0]+Shape[1]*b[1]+Shape[2]*b[2]+Shape[3]*b[3];
  N[2] = Shape[0]*c[0]+Shape[1]*c[1]+Shape[2]*c[2]+Shape[3]*c[3];

  return(0.25*sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

void
FaceQuad4::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape functions' derivatives w.r.t. the local coordinates
  double dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  dMdx[0] = 0.0; dMdx[1] = 0.0; dMdx[2] = 0.0;
  dMdy[0] = 0.0; dMdy[1] = 0.0; dMdy[2] = 0.0;

  dMdx[0] += dShapex[0]*X[0] + dShapex[1]*X[1] + dShapex[2]*X[2] + dShapex[3]*X[3];
  dMdx[1] += dShapex[0]*Y[0] + dShapex[1]*Y[1] + dShapex[2]*Y[2] + dShapex[3]*Y[3];
  dMdx[2] += dShapex[0]*Z[0] + dShapex[1]*Z[1] + dShapex[2]*Z[2] + dShapex[3]*Z[3];

  dMdy[0] += dShapey[0]*X[0] + dShapey[1]*X[1] + dShapey[2]*X[2] + dShapey[3]*X[3];
  dMdy[1] += dShapey[0]*Y[0] + dShapey[1]*Y[1] + dShapey[2]*Y[2] + dShapey[3]*Y[3];
  dMdy[2] += dShapey[0]*Z[0] + dShapey[1]*Z[1] + dShapey[2]*Z[2] + dShapey[3]*Z[3];
}

/*
void
FaceQuad4::IsoParamInterpolation(double* V, double* m, double* NdVals, int nComps=1, int* NdMap=0)
{
  double Shape[4];
  GetShapeFct(Shape,m);

  double (*ndVals)[nCmps] = reinterpret_cast<double (*)[nCmps]>(NdVals);  
  if(NdMap)
    for(int i=0; i<nComps; i++)
      V[i] = Shape[0]*ndVals[NdMap[0]][i]+Shape[1]*ndVals[NdMap[1]][i]
            +Shape[2]*ndVals[NdMap[2]][i]+Shape[3]*ndVals[NdMap[3]][i];
  else
    for(int i=0; i<nComps; i++)
      V[i] = Shape[0]*ndVals[0][i]+Shape[1]*ndVals[1][i]
            +Shape[2]*ndVals[2][i]+Shape[3]*ndVals[3][i];
}

void
IsoParamInterpolation(double* V, double* Shape, double* NdVals, int nNds, int nComps=1, int* NdMap=0)
{
  double (*ndVals)[nCmps] = reinterpret_cast<double (*)[nCmps]>(NdVals);  
  std::fill(V,V+nComps,0.0);
  if(NdMap)
    for(int j=0; j<nNds; ++j)
      for(int i=0; i<nComps; ++i) 
        V[i] += Shape[j]*ndVals[NdMap[j]][i];
  else
    for(int j=0; j<nNds; ++j)
      for(int i=0; i<nComps; ++i) 
        V[i] += Shape[j]*ndVals[j][i];
}
*/

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad4::LocalToGlobalCoord(double* M, double* m, CoordSet &cs)
{
  // input  : m, the local coordinates of some point P on the face
  //         cs, the global coordinates of the vertices of the face 
  // output : M, the global x,y,z coordinates of P
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);

  double Shape[4];
  GetShapeFct(Shape,m);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  M[0] = Shape[0]*X[0]+Shape[1]*X[1]+Shape[2]*X[2]+Shape[3]*X[3];
  M[1] = Shape[0]*Y[0]+Shape[1]*Y[1]+Shape[2]*Y[2]+Shape[3]*Y[3];
  M[2] = Shape[0]*Z[0]+Shape[1]*Z[1]+Shape[2]*Z[2]+Shape[3]*Z[3];
}

void 
FaceQuad4::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

#ifndef NICER_IMPLEMENTATION
double
FaceQuad4::GetJacobian(double *m, CoordSet &cs)
{
  double x = m[0];
  double y = m[1];

  double d1 = 0.5*(1.0+x);
  double d2 = 0.5*(1.0+y);
  double d3 = 1.0-d1;
  double d4 = 1.0-d2;

  double Shape[4];
  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;

  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  double a[4], b[4], c[4];
  a[0] = (Y[1]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[1]-Z[0]);
  a[1] = (Y[1]-Y[0])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[1]-Z[0]);
  a[2] = (Y[2]-Y[3])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[2]-Z[3]);
  a[3] = (Y[2]-Y[3])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[3]);

  b[0] = (Z[1]-Z[0])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[1]-X[0]);
  b[1] = (Z[1]-Z[0])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[1]-X[0]);
  b[2] = (Z[2]-Z[3])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[2]-X[3]);
  b[3] = (Z[2]-Z[3])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[2]-X[3]);

  c[0] = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]);
  c[1] = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]);
  c[2] = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]);
  c[3] = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]);

  double N[3];
  N[0] = Shape[0]*a[0]+Shape[1]*a[1]+Shape[2]*a[2]+Shape[3]*a[3];
  N[1] = Shape[0]*b[0]+Shape[1]*b[1]+Shape[2]*b[2]+Shape[3]*b[3];
  N[2] = Shape[0]*c[0]+Shape[1]*c[1]+Shape[2]*c[2]+Shape[3]*c[3];

  return(0.25*sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

double
FaceQuad4::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
{
  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);
 
  double Shape[4];
  GetShapeFct(Shape,m);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  double a[4], b[4], c[4];
  a[0] = (Y[1]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[1]-Z[0]); 
  a[1] = (Y[1]-Y[0])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[1]-Z[0]); 
  a[2] = (Y[2]-Y[3])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[2]-Z[3]); 
  a[3] = (Y[2]-Y[3])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[3]); 

  b[0] = (Z[1]-Z[0])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[1]-X[0]); 
  b[1] = (Z[1]-Z[0])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[1]-X[0]); 
  b[2] = (Z[2]-Z[3])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[2]-X[3]); 
  b[3] = (Z[2]-Z[3])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[2]-X[3]); 

  c[0] = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]); 
  c[1] = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]); 
  c[2] = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]); 
  c[3] = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]); 

  Normal[0] = Shape[0]*a[0]+Shape[1]*a[1]+Shape[2]*a[2]+Shape[3]*a[3];
  Normal[1] = Shape[0]*b[0]+Shape[1]*b[1]+Shape[2]*b[2]+Shape[3]*b[3];
  Normal[2] = Shape[0]*c[0]+Shape[1]*c[1]+Shape[2]*c[2]+Shape[3]*c[3];
  
  double NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(NormN!=0.0){
    Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
  }
  //cerr << "here in FaceQuad4::GetIsoParamMappingNormalAndJacobian, Normal = " << Normal[0] << " " << Normal[1] << " " << Normal[2] << ", J = " << 0.25*NormN << endl;
  return(0.25*NormN); 
}

#else // NICER IMPLEMENTATION ...

double
FaceQuad4::GetJacobian(double* m, CoordSet& cs)
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
FaceQuad4::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // N = dM/dx x dM/dy
  Normal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  Normal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  Normal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  // PJSA normalizing the Normal seems to be unnecessary since the normal is re-multiplied by it's norm later
  //      this just makes the derivatives more complicated to compute
  return 1.0;
/*
  double NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(NormN!=0.0){
   Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
  }
  //cerr << "here in FaceQuad4::GetIsoParamMappingNormalAndJacobian, Normal = " << Normal[0] << " " << Normal[1] << " " << Normal[2] << ", J = " << NormN << endl;
  return(NormN);
*/
}
#endif

void
FaceQuad4::GetdNormal(double dNormal[][3], double* m, CoordSet& cs)
{
  // This function computes ddNormal which is the Jacobian (matrix) of the Normal multiplied by the verticies' coordinates
  // It is used to compute the gradient of the gap function

  // Compute shape functions' derivatives w.r.t. the local coordinates
  double dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);
  Node &nd4 = cs.getNode(Nodes[3]);

  double X[4], Y[4], Z[4];
  X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
  X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
  X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;
  X[3] = nd4.x; Y[3] = nd4.y; Z[3] = nd4.z;

  double dNdX[3][4], dNdY[3][4], dNdZ[3][4];
  for(int i = 0; i < 4; ++ i) {
    dNormal[3*i  ][0] = 0;
    dNormal[3*i  ][1] = dMdx[2]*dShapey[i] - dShapex[i]*dMdy[2];
    dNormal[3*i  ][2] = dShapex[i]*dMdy[1] - dMdx[1]*dShapey[i];
    dNormal[3*i+1][0] = dShapex[i]*dMdy[2] - dMdx[2]*dShapey[i];
    dNormal[3*i+1][1] = 0;
    dNormal[3*i+1][2] = dMdx[0]*dShapey[i] - dShapex[i]*dMdy[0];
    dNormal[3*i+2][0] = dMdx[1]*dShapey[i] - dShapex[i]*dMdy[1];
    dNormal[3*i+2][1] = dShapex[i]*dMdy[0] - dMdx[0]*dShapey[i];
    dNormal[3*i+2][2] = 0;
  }
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceQuad4::ComputeArea(CoordSet &cs,const int ngp=2)
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
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceQuad4::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[4];
  double xi, eta, weight, m[2];

  FullM Mass(4);
  Mass.zero();

   for(int igp=0;igp<ngp;igp++){
     for(int jgp=0;jgp<ngp;jgp++){
       _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
       m[0] = xi; m[1] = eta;
       dA = GetShapeFctAndJacobian(Shape, m, cs);
       // upper part
       for(int i=0;i<4;i++)
	 for(int j=i;j<4;j++) 
           Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];     
    }
  }
  // lower part by symmetry 
  for(int i=0;i<4;i++)
     for(int j=0;j<i;j++) 
        Mass[i][j] = Mass[j][i];     

  return(Mass);
}

void 
FaceQuad4::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<4; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[4];
  double xi, eta, weight, m[2];

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
       m[0] = xi; m[1] = eta;
       dA = GetShapeFctAndJacobian(Shape, m, cs);
       for(int i=0;i<4;i++)
         ShapeIntg[i] += rho*weight*dA*Shape[i];
    }
  }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad4::printNodes()
{
  filePrint(stderr,"   # Quad4 face el., nodes = %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],Nodes[3]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad4::print()
{
  printNodes();
}
