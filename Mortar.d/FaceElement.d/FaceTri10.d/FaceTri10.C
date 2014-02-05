// ----------------------------------------------------------------
// HB - 05/05/05
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING 
//    s
//    ^
//    |
//   2+
//    |\             Shape functions
//    | \            ---------------
//   6+  +5            see: FaceTri10::GetShapeFct
//    |   \
//    | +  \
//   7+ 10  +4
//    |      \
//   3+-+---+-+1 -> r
//      8   9
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
#include <Mortar.d/FaceElement.d/FaceTri10.d/FaceTri10.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri10::RefCoords[10][2] = {{ 1.0  , 0.0   },
                                      { 0.0  , 1.0   },
                                      { 0.0  , 0.0   },
                                      { 2./3., 1./3. },
                                      { 1./3., 2./3. },
                                      { 0.0  , 2./3. },
                                      { 0.0  , 1./3. }, 
                                      { 1./3., 0.0   },
                                      { 2./3., 0.0   },
                                      { 1./3., 1./3. }};
double* 
FaceTri10::ViewRefCoords() { return(FaceTri10::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri10::FaceTri10(int* nodenums)
{
  for(int i=0; i<10; i++) { Nodes[i] = nodenums[i]; }
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  for(int i=0; i<10; i++) { Nodes[i] = OldToNewNodeIds[Nodes[i]]; }
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
int
FaceTri10::nTri3Nodes() { return 3; }

int
FaceTri10::GetTri3Node(int i) { return Nodes[i]; }

void
FaceTri10::GetTri3Nodes(int *p, int* renumTable)
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
FaceTri10::GetTri3Nodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri10::nNodes() { return(10); }

int
FaceTri10::GetNode(int i) { return(Nodes[i]); }

void
FaceTri10::GetNodes(int *p, int* renumTable)
{
  if(renumTable) {
    for(int i=0; i<10; i++) 
       p[i] = renumTable[Nodes[i]];
  } else {
    for(int i=0; i<10; i++) 
       p[i] = Nodes[i];
  }
}

void
FaceTri10::GetNodes(int *p, std::map<int,int>& renumTable)
{
  for(int i=0; i<10; i++) 
    p[i] = renumTable[Nodes[i]];
}

int 
FaceTri10::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<10; i++)
    if(gNode==Nodes[i]){ found = true; break; }
  if(!found)  
    filePrint(stderr," ### In FaceTri10::GetNodeIndex(): node () does not belong to this element\n",gNode);
  return(i); 
}

int
FaceTri10::nVertices() { return(nTri3Nodes()); }

int
FaceTri10::GetVertex(int i) { return(GetTri3Node(i)); }

void
FaceTri10::GetVertices(int *p, int* renumTable) { GetTri3Nodes(p, renumTable); }

void
FaceTri10::GetVertices(int *p, std::map<int,int>& renumTable) { GetTri3Nodes(p, renumTable); }

int
FaceTri10::GetFaceElemType() { return FaceElement::TRIFACEC10; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri10::GetACMEFaceElemType() { return ContactSearch::TRIFACEQ6; }
#else
int
FaceTri10::GetACMEFaceElemType() { return 4; }
#endif

// As ACME doesn't support the Tri10 face element for
// FaceFaceInteraction (), we will pass to it the Tri3 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Tri10 face element has straight edges, but its is an 
// APPROXIMATION in the general case (i.e. curved edges/face). 
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri10::GetACMEFFIFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int 
FaceTri10::GetACMEFFIFaceElemType() { return 3; }
#endif
// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void 
FaceTri10::GetShapeFct(double *Shape, double *m)
{
  double r = m[0];
  double s = m[1];
  double t = 1.-r-s; 
  double ninehalf= 9./2.;
  double twonine = 2./9.;

  // following ACME Tri3 ordering ...
  Shape[0] = ninehalf*r*(r*(r-1.)+twonine);
  Shape[1] = ninehalf*s*(s*(s-1.)+twonine);
  Shape[2] = ninehalf*t*(t*(t-1.)+twonine);

  Shape[3] = ninehalf*s*r*(3.*r-1.); 
  Shape[4] = ninehalf*s*r*(3.*s-1.);

  Shape[5] = ninehalf*t*s*(3.*s-1.); 
  Shape[6] = ninehalf*t*s*(3.*t-1.);

  Shape[7] = ninehalf*t*r*(3.*t-1.);
  Shape[8] = ninehalf*t*r*(3.*r-1.); 
  
  Shape[9] = 27.*r*s*t;
}

void
FaceTri10::GetdShapeFct(double *dShapex, double *dShapey, double *m)
{
  double r = m[0];
  double s = m[1];
  double t = 1.-r-s;
  double ninehalf= 9./2.;
  double twonine = 2./9.;

  // !! idem ACME !!
  dShapex[0] =  ninehalf*(r*(3.*r-2.)+twonine);
  dShapex[1] =  0.0;
  dShapex[2] = -ninehalf*(t*(3.*t-2.)+twonine);
  dShapex[3] =  ninehalf*s*(6.*r-1.);
  dShapex[4] =  ninehalf*s*(3.*s-1.);
  dShapex[5] = -ninehalf*s*(3.*s-1.);
  dShapex[6] = -ninehalf*s*(6.*t-1.);
  dShapex[7] =  ninehalf*( t*(3.*t-1.) - r*(6.*t-1.));
  dShapex[8] =  ninehalf*(-r*(3.*r-1.) + t*(6.*r-1.));
  dShapex[9] =  27.*s*(t - r);

  dShapey[0] =  0.0      ;
  dShapey[1] =  ninehalf*(s*(3.*s-2.)+twonine);
  dShapey[2] = -ninehalf*(t*(3.*t-2.)+twonine);
  dShapey[3] =  ninehalf*r*(3.*r-1.);
  dShapey[4] =  ninehalf*r*(6.*s-1.);
  dShapey[5] =  ninehalf*(-s*(3.*s-1.) + t*(6.*s-1.));
  dShapey[6] =  ninehalf*( t*(3.*t-1.) - s*(6.*t-1.));
  dShapey[7] = -ninehalf*r*(6.*t-1.);
  dShapey[8] = -ninehalf*r*(3.*r-1.);
  dShapey[9] =  27.*r*(t - s);
}

double
FaceTri10::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  GetShapeFct(Shape, m);  
  return(GetJacobian(m, cs));
}

void
FaceTri10::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape fcts derivatives
  double dShapex[10], dShapey[10];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  double X[10], Y[10], Z[10];
  for(int i=0; i<10; i+=2){
    Node &ndi = cs.getNode(Nodes[i]);
    X[i]  = ndi.x ; Y[i]  = ndi.y ; Z[i]  = ndi.z ;
    Node &ndii= cs.getNode(Nodes[i+1]);
    X[i+1]= ndii.x; Y[i+1]= ndii.y; Z[i+1]= ndii.z;
  }

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<10; i+=5){
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2] + dShapex[i+3]*X[i+3] + dShapex[i+4]*X[i+4];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2] + dShapex[i+3]*Y[i+3] + dShapex[i+4]*Y[i+4];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2] + dShapex[i+3]*Z[i+3] + dShapex[i+4]*Z[i+4];
    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2] + dShapey[i+3]*X[i+3] + dShapey[i+4]*X[i+4];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2] + dShapey[i+3]*Y[i+3] + dShapey[i+4]*Y[i+4];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2] + dShapey[i+3]*Z[i+3] + dShapey[i+4]*Z[i+4];
  }
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  double Shape[10];
  GetShapeFct(Shape,m);

  double X[10], Y[10], Z[10];
  for(int i=0; i<10; i+=2){
    Node &ndi = cs.getNode(Nodes[i]);
    X[i]  = ndi.x ; Y[i]  = ndi.y ; Z[i]  = ndi.z ;
    Node &ndii= cs.getNode(Nodes[i+1]);
    X[i]  = ndi.x ; Y[i]  = ndi.y ; Z[i]  = ndi.z ;
  }

  M[0] = M[1] = M[2] = 0.0;
  for(int i=0; i<10; i+=5){
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2] + Shape[i+3]*X[i+3] + Shape[i+4]*X[i+4];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2] + Shape[i+3]*Y[i+3] + Shape[i+4]*Y[i+4];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2] + Shape[i+3]*Z[i+3] + Shape[i+4]*Z[i+4];
  }
}

void 
FaceTri10::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

double
FaceTri10::GetJacobian(double *m, CoordSet &cs)
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
FaceTri10::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
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
FaceTri10::GetIsoParamMappingNormalJacobianProduct(double *JNormal, double *m, CoordSet &cs)
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
FaceTri10::ComputeArea(CoordSet &cs,const int ngp=0)
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
FaceTri10::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass;
  Mass.zero();
  fprintf(stderr," *** WARNING: FaceTri10::ScalarMass(): NOT IMPLEMENTED. Return zero mass matrix.\n");
  return(Mass);
}

void 
FaceTri10::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  fprintf(stderr," *** WARNING: FaceTri10::IntegrateShapeFcts(): NOT IMPLEMENTED. Return zero shape functions'integral.\n");
  for(int i=0;i<10;i++) { ShapeIntg[i] = 0.0; }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceTri10::printNodes()
{
  filePrint(stderr,"   # Tri10 face el., nodes = %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
           Nodes[0],Nodes[1],Nodes[2],Nodes[3],Nodes[4],Nodes[5],Nodes[6],Nodes[7],Nodes[8],Nodes[9]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::print()
{
  printNodes();
}
