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
#include <cstdio>
#include <cstdlib>
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
  return LocalToGlobalCoordTemp(M, m, cs);
}

#if (MAX_MORTAR_DERIVATIVES > 0)
void
FaceTri3::LocalToGlobalCoord(ActiveDouble *M, ActiveDouble *m, MadCoordSet &cs)
{
  return LocalToGlobalCoordTemp(M, m, cs);
}
#endif

void 
FaceTri3::GetShapeFctVal(double *Shape, double *m)
{
   GetShapeFct(Shape, m);
}

#if (MAX_MORTAR_DERIVATIVES > 0)
void
FaceTri3::GetShapeFctVal(ActiveDouble *Shape, ActiveDouble *m)
{
   GetShapeFct(Shape, m);
}
#endif

double
FaceTri3::GetJacobian(double *m, CoordSet &cs)
{
   return(GetJacobian<double,CoordSet>(cs));
}

double
FaceTri3::GetIsoParamMappingNormalAndJacobian(double *Normal, double *m, CoordSet &cs)
{
   GetIsoParamMappingNormalJacobianProduct(Normal, m, cs);

   double NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
   
   if(NormN!=0.0){
     Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
   }
   return(NormN); // !! A CONTROLER !! 
}

void
FaceTri3::GetIsoParamMappingNormalJacobianProduct(double *JNormal, double *m, CoordSet &cs)
{
   GetIsoParamMappingNormalJacobianProductTemp(JNormal, m, cs);
}

#if (MAX_MORTAR_DERIVATIVES > 0)
void
FaceTri3::GetIsoParamMappingNormalJacobianProduct(ActiveDouble *JNormal, ActiveDouble *m, MadCoordSet &cs)
{
   GetIsoParamMappingNormalJacobianProductTemp(JNormal, m, cs);
}
#endif

void
FaceTri3::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  // This function computes dJNormal which is the Jacobian (matrix) of J times the unit normal 
  // It is used to compute the gradient of the gap function

  // Compute shape functions' derivatives w.r.t. the local coordinates
  double dShapex[3], dShapey[3];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  double dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute dJNormal
  for(int i = 0; i < 3; ++ i) {
    dJNormal[3*i  ][0] = 0;
    dJNormal[3*i  ][1] = dMdx[2]*dShapey[i] - dShapex[i]*dMdy[2];
    dJNormal[3*i  ][2] = dShapex[i]*dMdy[1] - dMdx[1]*dShapey[i];
    dJNormal[3*i+1][0] = dShapex[i]*dMdy[2] - dMdx[2]*dShapey[i];
    dJNormal[3*i+1][1] = 0;
    dJNormal[3*i+1][2] = dMdx[0]*dShapey[i] - dShapex[i]*dMdy[0];
    dJNormal[3*i+2][0] = dMdx[1]*dShapey[i] - dShapex[i]*dMdy[1];
    dJNormal[3*i+2][1] = dShapex[i]*dMdy[0] - dMdx[0]*dShapey[i];
    dJNormal[3*i+2][2] = 0;
  }
/*

//  std::cerr << "dNormal[][0] = "; for(int i=0; i<9; ++i) std::cerr << dNormal[i][0] << " "; std::cerr << std::endl;
//  std::cerr << "dNormal[][1] = "; for(int i=0; i<9; ++i) std::cerr << dNormal[i][1] << " "; std::cerr << std::endl;
//  std::cerr << "dNormal[][2] = "; for(int i=0; i<9; ++i) std::cerr << dNormal[i][2] << " "; std::cerr << std::endl;

  Node &nd1 = cs.getNode(Nodes[0]);
  Node &nd2 = cs.getNode(Nodes[1]);
  Node &nd3 = cs.getNode(Nodes[2]);

  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  x1 = nd1.x; y1 = nd1.y; z1 = nd1.z;
  x2 = nd2.x; y2 = nd2.y; z2 = nd2.z;
  x3 = nd3.x; y3 = nd3.y; z3 = nd3.z;

//  This is an alternative expression for dJNormal
//  double dN0[9][3] = {
//    {       0, z3 - z2, y2 - y3},
//    { z2 - z3,       0, x3 - x2},
//    { y3 - y2, x2 - x3,       0},
//    {       0, z1 - z3, y3 - y1},
//    { z3 - z1,       0, x1 - x3},
//    { y1 - y3, x3 - x1,       0},
//    {       0, z2 - z1, y1 - y2},
//    { z1 - z2,       0, x2 - x1},
//    { y2 - y1, x1 - x2,       0} };

//  std::cerr << "dN0[][0] =     "; for(int i=0; i<9; ++i) std::cerr << dN0[i][0] << " "; std::cerr << std::endl;
//  std::cerr << "dN0[][1] =     "; for(int i=0; i<9; ++i) std::cerr << dN0[i][1] << " "; std::cerr << std::endl;
//  std::cerr << "dN0[][2] =     "; for(int i=0; i<9; ++i) std::cerr << dN0[i][2] << " "; std::cerr << std::endl;

  for(int i=0;i<9;++i)
    for(int j=0;j<3;++j)
      dJNormal[i][j] = dN0[i][j];
*/
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
   double Area = 0.5*GetJacobian<double,CoordSet>(cs);
   Area *= rho/24;
   Mass[0][0] = 2.*Area; Mass[0][1] =    Area; Mass[0][1] =    Area;
   Mass[1][0] =    Area; Mass[1][1] = 2.*Area; Mass[1][1] =    Area;
   Mass[2][0] =    Area; Mass[2][1] =    Area; Mass[2][2] = 2.*Area;

   return(Mass);
}

void 
FaceTri3::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
   double Area = 0.5*GetJacobian<double,CoordSet>(cs);
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
int* FaceTri3::dofs(DofSetArray &dsa, int *p, int *fnId) 
{
  if(p == 0) p = new int[9];
    dsa.number(fnId[Nodes[0]], DofSet::XYZdisp, p);
    dsa.number(fnId[Nodes[1]], DofSet::XYZdisp, p+3);
    dsa.number(fnId[Nodes[2]], DofSet::XYZdisp, p+6);
    return p;
}

void FaceTri3::computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, 
                           GeomState*, int *fnId) 
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(fnId[Nodes[0]], xyz[0], xyz[0]+3);
  state.getDV(fnId[Nodes[1]], xyz[1], xyz[1]+3);
  state.getDV(fnId[Nodes[2]], xyz[2], xyz[2]+3);

  for(int j=0; j<6; ++j)
    res[j] = gp[0]*xyz[0][j] + gp[1]*xyz[1][j] + (1.0-gp[0]-gp[1])*xyz[2][j]; //using ACME convention
}

void FaceTri3::getFlLoad(const InterpPoint &ip, double *flF, double *resF) 
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]    = gp[0] * flF[i]; //using ACME convention
    resF[3+i]  = gp[1] * flF[i];
    resF[6+i] = (1.0-gp[0]-gp[1]) * flF[i];
  }
}
