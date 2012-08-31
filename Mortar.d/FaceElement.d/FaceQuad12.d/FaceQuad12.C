// ----------------------------------------------------------------
// HB - 05/10/05
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//          y
//          ^
//          |           Shape functions
//       10   9         ---------------
//   4+---+---+---+3      Phi[1] = (1/4).(1-x).(1-y).(-1-x-y)
//    |           |       Phi[2] = (1/4).(1+x).(1-y).(-1+x-y) 
//  11+           +8      Phi[3] = (1/4).(1+x).(1+y).(-1+x+y) 
//    |           | -> x  Phi[4] = (1/4).(1-x).(1+y).(-1-x+y) 
//  12+           +7      Phi[5] = (1/2).(1-x.x).(1-y) 
//    |           |       Phi[6] = (1/2).(1+x).(1-y.y) 
//   1+---+---+---+2      Phi[7] = (1/2).(1-x.x).(1+y) 
//        5   6           Phi[8] = (1/2).(1-x).(1-y.y) 
//                        Phi[12] =       (1-x.x).(1-y.y) 
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
#include <Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.h>

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
double FaceQuad12::RefCoords[12][2] = {{-1.0  ,-1.0  },
                                       { 1.0  ,-1.0  },
                                       { 1.0  , 1.0  },
                                       {-1.0  , 1.0  },
                                       {-1./3.,-1.0  },
                                       { 1./3.,-1.0  },
                                       { 1.0  ,-1./3.},
                                       { 1.0  , 1./3.},
                                       { 1./3., 1.0  },
                                       {-1./3., 1.0  },
                                       {-1.0  , 1./3.},
                                       {-1.0  ,-1./3.}};
double* 
FaceQuad12::ViewRefCoords() { return(FaceQuad12::RefCoords[0]); }


// -----------------------------------------------------------------------------------------------------
//                              CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad12::FaceQuad12(int* nodenums)
{
  for(int i=0; i<12; i++) { Nodes[i] = nodenums[i]; }
}

// copy constructor
/*
FaceQuad12::FaceQuad12(const FaceQuad12 &FQ8)
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
FaceQuad12::clone()
{
  return new FaceQuad12(*this);
}
*/

void 
FaceQuad12::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  for(int i=0; i<12; i++) { Nodes[i] = OldToNewNodeIds[Nodes[i]]; }
}

// -----------------------------------------------------------------------------------------------------
//                              GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF LOCAL METHODS
// -------------------------------
int
FaceQuad12::nQuad4Nodes() { return 4; }

int
FaceQuad12::GetQuad4Node(int i) { return(Nodes[i]); } 

void
FaceQuad12::GetQuad4Nodes(int *p, int* renumTable)
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
FaceQuad12::GetQuad4Nodes(int *p, std::map<int,int>&renumTable)
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
FaceQuad12::nNodes() { return(12); }

int
FaceQuad12::GetNode(int i) { return(Nodes[i]); }

void
FaceQuad12::GetNodes(int *p, int* renumTable)
{
  if(renumTable) {
    for(int i=0; i<12; i++)
      p[i] = renumTable[Nodes[i]];
  } else {
    for(int i=0; i<12; i++)
      p[i] = Nodes[i];
  }
}

void
FaceQuad12::GetNodes(int *p, std::map<int,int>& renumTable)
{
  for(int i=0; i<12; i++)
    p[i] = renumTable[Nodes[i]];
}

int 
FaceQuad12::GetNodeIndex(int gNode)
{
  int i;
  bool found = false; 
  for(i=0; i<12; i++)
    if(gNode==Nodes[i]){ found = true; break; }
  if(!found)  
    filePrint(stderr," *** WARNING: FaceQuad12::GetNodeIndex(...): node (%6d) does not belong to this element\n",gNode);
  return(i); 
}

int
FaceQuad12::GetFaceElemType() { return(FaceElement::QUADFACEC12); }

#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad12::GetACMEFaceElemType() { return ContactSearch::QUADFACEQ8; }
#else
 int
 FaceQuad12::GetACMEFaceElemType() { return 2; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int 
FaceQuad12::nVertices() { return(nQuad4Nodes()); }

int 
FaceQuad12::GetVertex(int i) { return(GetQuad4Node(i)); }

void 
FaceQuad12::GetVertices(int* p, int* renumTable) { GetQuad4Nodes(p, renumTable); }

void 
FaceQuad12::GetVertices(int* p, std::map<int,int>&renumTable) { GetQuad4Nodes(p, renumTable); }

// As ACME doesn't support the Quad12 face element for
// FaceFaceInteraction (FFI), we will pass to it the Quad4 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Quad12 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
 ContactSearch::ContactFace_Type
 FaceQuad12::GetACMEFFIFaceElemType() { return(ContactSearch::QUADFACEL4); }
#else
 int
 FaceQuad12::GetACMEFFIFaceElemType() { return 1; }
#endif

// -----------------------------------------------------------------------------------------------------
//                          MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad12::GetShapeFct(double *Shape, double *m)
{
  double x = m[0];
  double y = m[1];
  double cv  = 1./32.;
  double cm  = 9./32.;

  double xm = 1.-x;
  double xp = 1.+x;
  double ym = 1.-y;
  double yp = 1.+y;
  double d2 = 10.-9.*(x*x+y*y);

  Shape[ 0] = -cv*xm*ym*d2;
  Shape[ 1] = -cv*xp*ym*d2;
  Shape[ 2] = -cv*xp*yp*d2;
  Shape[ 3] = -cv*xm*yp*d2;

  double dm3x = 1.-3.*x;
  double dp3x = 1.+3.*x;
  double dm3y = 1.-3.*y;
  double dp3y = 1.+3.*y;

  Shape[ 4] = cm*xm*xp*ym*dm3x;
  Shape[ 5] = cm*xm*xp*ym*dp3x;
  Shape[ 6] = cm*xp*ym*yp*dm3y;
  Shape[ 7] = cm*xp*ym*yp*dp3y;
  Shape[ 8] = cm*xm*xp*yp*dp3x;
  Shape[ 9] = cm*xm*xp*yp*dm3x;
  Shape[10] = cm*xm*ym*yp*dp3y;
  Shape[11] = cm*xm*ym*yp*dm3y;
}

void
FaceQuad12::GetdShapeFct(double *dShapex, double *dShapey, double *m)
{
  double x = m[0];
  double y = m[1];
  double cv  = 1./32.;
  double cm  = 9./32.;

  double xm  = 1.-x;
  double xp  = 1.+x;
  double ym  = 1.-y;
  double yp  = 1.+y;
  double d2  = 10.-9.*(x*x+y*y);
  double dd2x= 18.*x;
  double dd2y= 18.*y;

  dShapex[0] = cv*ym*( d2 + xm*dd2x);
  dShapex[1] = cv*ym*(-d2 + xp*dd2x);
  dShapex[2] = cv*yp*(-d2 + xp*dd2x);
  dShapex[3] = cv*yp*( d2 + xm*dd2x);

  dShapey[0] = cv*xm*( d2 + ym*dd2y);
  dShapey[1] = cv*xp*( d2 + ym*dd2y);
  dShapey[2] = cv*xp*(-d2 + yp*dd2y);
  dShapey[3] = cv*xm*(-d2 + yp*dd2y);

  double dm3x = 1.-3.*x;
  double dp3x = 1.+3.*x;
  double dm3y = 1.-3.*y;
  double dp3y = 1.+3.*y;

  dShapex[ 4] =  cm*ym*(-3.-x*(2-9.*x));  
  dShapex[ 5] =  cm*ym*( 3.-x*(2+9.*x));  
  dShapex[ 6] =  cm*ym*yp*dm3y;  
  dShapex[ 7] =  cm*ym*yp*dp3y;  
  dShapex[ 8] =  cm*yp*( 3.-x*(2+9.*x));  
  dShapex[ 9] =  cm*yp*(-3.-x*(2-9.*x));  
  dShapex[10] = -cm*ym*yp*dp3y;  
  dShapex[11] = -cm*ym*yp*dm3y;  

  dShapey[ 4] = -cm*xm*xp*dm3x;  
  dShapey[ 5] = -cm*xm*xp*dp3x;  
  dShapey[ 6] =  cm*xp*(-3.-y*(2-9.*y));  
  dShapey[ 7] =  cm*xp*( 3.-y*(2+9.*y));  
  dShapey[ 8] =  cm*xm*xp*dp3x;  
  dShapey[ 9] =  cm*xm*xp*dm3x;  
  dShapey[10] =  cm*xm*( 3.-y*(2+9.*y));  
  dShapey[11] =  cm*xm*(-3.-y*(2-9.*y));  
}

double
FaceQuad12::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
  GetShapeFct(Shape, m);
  return(GetJacobian(m, cs));
}

void
FaceQuad12::ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs)
{
  // Compute shape fcts derivatives
  double dShapex[12], dShapey[12];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  double X[12], Y[12], Z[12];
  for(int i=0; i<12; i+=2){
    Node &ndi = cs.getNode(Nodes[i]);
    X[i]  = ndi.x ; Y[i]  = ndi.y ; Z[i]  = ndi.z ;
    Node &ndii= cs.getNode(Nodes[i+1]);
    X[i+1]= ndii.x; Y[i+1]= ndii.y; Z[i+1]= ndii.z;
  }

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<12; i+=4){
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2] + dShapex[i+3]*X[i+3];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2] + dShapex[i+3]*Y[i+3];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2] + dShapex[i+3]*Z[i+3];
    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2] + dShapey[i+3]*X[i+3];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2] + dShapey[i+3]*Y[i+3];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2] + dShapey[i+3]*Z[i+3];
  }
}

double
FaceQuad12::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
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

void
FaceQuad12::GetIsoParamMappingNormalJacobianProduct(double *JNormal, double *m, CoordSet &cs)
{
  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // JN = dM/dx x dM/dy
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad12::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  double Shape[12];
  GetShapeFct(Shape,m);

  double X[12], Y[12], Z[12];
  for(int i=0; i<12; i+=2){
    Node &ndi = cs.getNode(Nodes[i]);
    X[i]  = ndi.x ; Y[i]  = ndi.y ; Z[i]  = ndi.z ;
    Node &ndii= cs.getNode(Nodes[i+1]);
    X[i+1]= ndii.x; Y[i+1]= ndii.y; Z[i+1]= ndii.z;
  }

  M[0] = M[1] = M[2] = 0.0;
  for(int i=0; i<12; i+=4){
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2] + Shape[i+3]*X[i+3];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2] + Shape[i+3]*Y[i+3];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2] + Shape[i+3]*Z[i+3];
  }
}

void 
FaceQuad12::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFct(Shape, m);
}

double
FaceQuad12::GetJacobian(double *m, CoordSet &cs)
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
FaceQuad12::ComputeArea(CoordSet &cs,const int ngp=2)
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
FaceQuad12::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[12];
  double weight, m[2];

  FullM Mass(12);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,m[0],m[1],weight);
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<12;i++)
        for(int j=i;j<12;j++) 
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];    
    }
  }
  // lower part by symmetry 
  for(int i=0;i<12;i++)
    for(int j=0;j<i;j++) 
      Mass[i][j] = Mass[j][i];    

  return(Mass);
}

void 
FaceQuad12::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<12; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[12];
  double weight, m[2];

  for(int igp=0;igp<ngp;igp++){
    for(int jgp=0;jgp<ngp;jgp++){
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,m[0],m[1],weight);
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<12;i++)
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
FaceQuad12::printNodes()
{
  filePrint(stderr,"   # Quad12 face el., nodes = ");
  for(int i=0; i<12; i++)
    filePrint(stderr," %6d",Nodes[i]);
  filePrint(stderr,"\n");
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad12::print()
{
  printNodes();
}
