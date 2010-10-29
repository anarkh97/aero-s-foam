/* --------------------------------------------------------------
 HB - 08/15/03
 ----------------------------------------------------------------
 WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
    s 
    ^
    | 
   2+
    |\             Shape functions
    | \            ---------------
    |  \             Phi[1] = r
    |   \            Phi[2] = s
    |    \           Phi[3] = t = 1.-r-s
    |     \
   3+------+1 -> r

 ----------------------------------------------------------------*/

// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// STL
#include <map>

// FEM headers 
#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>
#include <Utils.d/dofset.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/State.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri3::RefCoords[3][2] = {{ 1.0 , 0.0 },
                                    { 0.0 , 1.0 },
                                    { 0.0 , 0.0 }};
double* 
FaceTri3::ViewRefCoords() { return(FaceTri3::RefCoords[0]); }

/*
double*
FaceTri3::ViewRefCoord(int i) { return(&(RefCoords[i][0])); }
*/
// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri3::FaceTri3(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void 
FaceTri3::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri3::nNodes() { return 3; }

int
FaceTri3::GetNode(int i) { return Nodes[i]; }

void
FaceTri3::GetNodes(int *p, int* renumTable)
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
FaceTri3::GetNodes(int *p, std::map<int,int>& renumTable)
{
   p[0] = renumTable[Nodes[0]];
   p[1] = renumTable[Nodes[1]];
   p[2] = renumTable[Nodes[2]];
}

int 
FaceTri3::GetNodeIndex(int gNode)
{
   int i;
   bool found = false; 
   for(i=0; i<3; i++)
     if(gNode==Nodes[i]){ found = true; break; }
   if(!found){  
     filePrint(stderr," *** WARNING: FaceTri3::GetNodeIndex(): node (%6d) does not belong to this element\n",gNode); 
     printNodes();
   }
   return(i); 
}

int
FaceTri3::GetFaceElemType() { return FaceElement::TRIFACEL3; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri3::GetACMEFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int
FaceTri3::GetACMEFaceElemType() { return 3; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceTri3::nVertices() { return nNodes(); }

int
FaceTri3::GetVertex(int i) { return GetNode(i); }

void
FaceTri3::GetVertices(int* p, int* renumTable) { GetNodes(p, renumTable); }

void
FaceTri3::GetVertices(int* p, std::map<int,int>& renumTable) { GetNodes(p, renumTable); }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri3::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#else
int
FaceTri3::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#endif
// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------

void
FaceTri3::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
   Node &nd1 = cs.getNode(Nodes[0]);
   Node &nd2 = cs.getNode(Nodes[1]);
   Node &nd3 = cs.getNode(Nodes[2]);

   double r = m[0];
   double s = m[1];
   double t = 1.-r-s; 

   double X[3], Y[3], Z[3];
   X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
   X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
   X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;

   //M[0] = t*X[0]+r*X[1]+s*X[2];
   //M[1] = t*Y[0]+r*Y[1]+s*Y[2];
   //M[2] = t*Z[0]+r*Z[1]+s*Z[2];
   // !! idem ACME !!
   M[0] = r*X[0]+s*X[1]+t*X[2];
   M[1] = r*Y[0]+s*Y[1]+t*Y[2];
   M[2] = r*Z[0]+s*Z[1]+t*Z[2];
}

void 
FaceTri3::GetShapeFct(double *Shape, double *m)
{
   double r = m[0];
   double s = m[1];
   double t = 1.-r-s; 

   //Shape[0] = t;
   //Shape[1] = r;
   //Shape[2] = s;
   // !! idem ACME !! 
   Shape[0] = r;
   Shape[1] = s;
   Shape[2] = t;
}

void 
FaceTri3::GetShapeFctVal(double *Shape, double *m)
{
   GetShapeFct(Shape, m);
}

double
FaceTri3::GetShapeFctAndJacobian(double *Shape, double *m, CoordSet &cs)
{
   double r = m[0];
   double s = m[1];
   double t = 1.-r-s; 

   //Shape[0] = t;
   //Shape[1] = r;
   //Shape[2] = s;
   // !! idem ACME !! 
   Shape[0] = r;
   Shape[1] = s;
   Shape[2] = t;

   return(GetJacobian(cs));
}

double
FaceTri3::GetJacobian(double *m, CoordSet &cs)
{
   return(GetJacobian(cs));
}

double
FaceTri3::GetJacobian(CoordSet &cs)
{
   // J = 2*Area = ||12 x 13||
   Node &nd1 = cs.getNode(Nodes[0]);
   Node &nd2 = cs.getNode(Nodes[1]);
   Node &nd3 = cs.getNode(Nodes[2]);

   double X[3], Y[3], Z[3];
   X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
   X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
   X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;

   double V12[3], V13[3];
   V12[0] = X[1]-X[0]; V12[1] = Y[1]-Y[0]; V12[2] = Z[1]-Z[0];
   V13[0] = X[2]-X[0]; V13[1] = Y[2]-Y[0]; V13[2] = Z[2]-Z[0];

   double N[3];
   N[0] = V12[1]*V13[2] - V12[2]*V13[1];
   N[1] = V12[2]*V13[0] - V12[0]*V13[2];
   N[2] = V12[0]*V13[1] - V12[1]*V13[0];

   return(sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

double
FaceTri3::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
{
   // J = 2*Area = ||12 x 13||
   Node &nd1 = cs.getNode(Nodes[0]);
   Node &nd2 = cs.getNode(Nodes[1]);
   Node &nd3 = cs.getNode(Nodes[2]);

   double X[3], Y[3], Z[3];
   X[0] = nd1.x; Y[0] = nd1.y; Z[0] = nd1.z;
   X[1] = nd2.x; Y[1] = nd2.y; Z[1] = nd2.z;
   X[2] = nd3.x; Y[2] = nd3.y; Z[2] = nd3.z;

   double V12[3], V13[3];
   V12[0] = X[1]-X[0]; V12[1] = Y[1]-Y[0]; V12[2] = Z[1]-Z[0];
   V13[0] = X[2]-X[0]; V13[1] = Y[2]-Y[0]; V13[2] = Z[2]-Z[0];

   Normal[0] = V12[1]*V13[2] - V12[2]*V13[1];
   Normal[1] = V12[2]*V13[0] - V12[0]*V13[2];
   Normal[2] = V12[0]*V13[1] - V12[1]*V13[0];

   double NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
   
   if(NormN!=0.0){
      Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
   }
   return(NormN); // !! A CONTROLER !! 
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceTri3::ComputeArea(CoordSet &cs,const int ngp=0)
{
  return 0.5*GetJacobian(cs);
} 
*/

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
FullM
FaceTri3::ScalarMass(CoordSet &cs, double rho, int ngp)
{
   FullM Mass;
   Mass.zero();
   double Area = 0.5*GetJacobian(cs);
   Area *= rho/24;
   Mass[0][0] = 2.*Area; Mass[0][1] =    Area; Mass[0][1] =    Area;
   Mass[1][0] =    Area; Mass[1][1] = 2.*Area; Mass[1][1] =    Area;
   Mass[2][0] =    Area; Mass[2][1] =    Area; Mass[2][2] = 2.*Area;

   return(Mass);
}

void 
FaceTri3::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
   double Area = 0.5*GetJacobian(cs);
   Area *= rho/6;
   ShapeIntg[0] = Area;
   ShapeIntg[1] = Area;
   ShapeIntg[2] = Area;
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
void
FaceTri3::printNodes()
{
  filePrint(stderr,"   # Tri3  face el., nodes = %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2]);
}

void
FaceTri3::print()
{
  printNodes();
}
// -----------------------------------------------------------------------------------------------------
//                                            FS COMMUNICATION (KW) 
// -----------------------------------------------------------------------------------------------------
int* FaceTri3::dofs(DofSetArray &dsa, int *p) 
{
  if(p == 0) p = new int[9];
    dsa.number(Nodes[0], DofSet::XYZdisp, p);
    dsa.number(Nodes[1], DofSet::XYZdisp, p+3);
    dsa.number(Nodes[2], DofSet::XYZdisp, p+6);
    return p;
}

void FaceTri3::computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, GeomState*) 
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(Nodes[0], xyz[0], xyz[0]+3);
  state.getDV(Nodes[1], xyz[1], xyz[1]+3);
  state.getDV(Nodes[2], xyz[2], xyz[2]+3);

  for(int j=0; j<6; ++j)
    res[j] = gp[0]*xyz[0][j] + gp[1]*xyz[1][j] + (1.0-gp[0]-gp[1])*xyz[2][j]; //using ACME convention
}

void FaceTri3::getFlLoad(CoordSet&cs, const InterpPoint &ip, double *flF, double *resF, GeomState*) 
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]    = gp[0] * flF[i]; //using ACME convention
    resF[3+i]  = gp[1] * flF[i];
    resF[6+i] = (1.0-gp[0]-gp[1]) * flF[i];
  }
}








