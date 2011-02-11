// ----------------------------------------------------------------
// HB - 05/06/03
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
#include <Mortar.d/FaceElement.d/FacePoint1.d/FacePoint1.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FacePoint1::RefCoords[1][2] = {{0.0,0.0}};

double*
FacePoint1::ViewRefCoords() { return(FacePoint1::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FacePoint1::FacePoint1(int* nodenums)
{
  Nodes[0] = nodenums[0];
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE  METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void 
FacePoint1::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FacePoint1::nNodes() { return 1; }

int
FacePoint1::GetNode(int i) { return Nodes[i]; }

void
FacePoint1::GetNodes(int *p, int* renumTable)
{
  if(renumTable){
    p[0] = renumTable[Nodes[0]];
  } else {
    p[0] = Nodes[0];
  }
}

void
FacePoint1::GetNodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
}

int 
FacePoint1::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<1; i++)
    if(gNode==Nodes[i]){ found = true; break; }
   if(!found){  
    filePrint(stderr," *** WARNING: FacePoint1::GetNodeIndex(): node (%6d) does not belong to this element\n",gNode);
    printNodes();
   }
  return i; 
}

int
FacePoint1::GetFaceElemType() { return FaceElement::POINTFACE; }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FacePoint1::GetACMEFaceElemType() { cerr << "POINTFACE not supported\n"; return (ContactSearch::ContactFace_Type) 0; }
#else
 int
 FacePoint1::GetACMEFaceElemType() { cerr << "POINTFACE not supported\n"; return 0; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int 
FacePoint1::nVertices() { return nNodes(); }

int 
FacePoint1::GetVertex(int i) { return GetNode(i); }

void 
FacePoint1::GetVertices(int* p, int* renumTable) { GetNodes(p, renumTable); }

void 
FacePoint1::GetVertices(int* p, std::map<int,int>& renumTable) { GetNodes(p, renumTable); }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FacePoint1::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#else
 int
 FacePoint1::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FacePoint1::GetShapeFct(double *Shape, double *m)
{
  cerr << "POINTFACE not supported\n"; 
}

void
FacePoint1::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  cerr << "POINTFACE not supported\n"; 
}

double
FacePoint1::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  cerr << "POINTFACE not supported\n"; return 0.0;
}

void
FacePoint1::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  cerr << "POINTFACE not supported\n"; 
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FacePoint1::LocalToGlobalCoord(double* M, double* m, CoordSet &cs)
{
  cerr << "POINTFACE not supported\n";
}

void 
FacePoint1::GetShapeFctVal(double *Shape, double *m)
{
  cerr << "POINTFACE not supported\n";
}

double
FacePoint1::GetJacobian(double *m, CoordSet &cs)
{
  cerr << "POINTFACE not supported\n"; return 0.0;
}

double
FacePoint1::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
{
  cerr << "POINTFACE not supported\n"; return 0.0;
}

FullM
FacePoint1::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass(1);
  Mass.zero();
  cerr << "POINTFACE not supported\n"; return 0;
  return(Mass);
}

void 
FacePoint1::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  cerr << "POINTFACE not supported\n";
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FacePoint1::printNodes()
{
  filePrint(stderr,"   # Point1 face el., nodes = %6d\n",Nodes[0]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FacePoint1::print()
{
  printNodes();
}
