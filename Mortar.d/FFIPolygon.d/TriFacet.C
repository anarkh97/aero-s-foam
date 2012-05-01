// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <Utils.d/dbg_alloca.h>
#include <cmath>

// FEM headers
#include <Element.d/Element.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

// External routines
extern void getGaussPtOnTriangle(int, int, double&, double&, double&, double&);

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
TriFacet<Scalar>::TriFacet()
{
  Initialize();
}

template<class Scalar>
TriFacet<Scalar>::TriFacet(FaceElement* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3)
{
  // set area to zero
  //Area   = 0.0;
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];

  // set ptr to face el.
  FaceEl = FaceElem;
}

// -----------------------------------------------------------------------------------------------------
//                                             DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
TriFacet<Scalar>::~TriFacet()
{
  // set area to zero
  //Area   = 0.0;
  // set ptr to face el. to NULL
  FaceEl = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
TriFacet<Scalar>::Initialize()
{
  // set area to zero
  //Area   = 0.0;
  // set local coord. of triangular vertices in face el. to 0.
  LocalCoordOnFaceEl[0][0] = 0.0; LocalCoordOnFaceEl[0][1] = 0.0;
  LocalCoordOnFaceEl[1][0] = 0.0; LocalCoordOnFaceEl[1][1] = 0.0;
  LocalCoordOnFaceEl[2][0] = 0.0; LocalCoordOnFaceEl[2][1] = 0.0;
  // set ptr to face el. to NULL
  FaceEl = 0;
}
// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
TriFacet<Scalar>::SetTriFacet(FaceElement* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3)
{
  // set area to zero
  //Area = 0.0;
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];
  // set ptr to face el.
  FaceEl = FaceElem;
}

//void
//TriFacet::SetArea(double S) { Area = S; };

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
TriFacet<Scalar>::Print()
{
  //fprintf(stderr," * TriFacet area = %e\n", Area);
  std::cerr << " * TriFacet vertices coord. (on ref. face elem.):\n";
  for(int i=0; i<3; i++)
    std::cerr << "  -> vertex " << i+1 << ": x = " << LocalCoordOnFaceEl[i][0]
              << ", y = " << LocalCoordOnFaceEl[i][1] << std::endl;
  std::cerr << " * TriFacet ptr to face elem " << FaceEl << std::endl;
}

/*
void
TriFacet::PrintVerticesXYZ(FILE* file=stderr, CoordSet& cs, int& firstVertId)
{
  double mOnFaceEl[2];
  double XYZ[3];
  for(int i=0; i<3; i++){
     double m[2] = {LocalCoordOnFaceEl[i][0],LocalCoordOnFaceEl[i][0]};
     LocalToLocalCoordOnFaceEl(m,mOnFaceEl);    
     FaceEl->LocalToGlobalCoord(XYZ,mOnFaceEl,cs);
     fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
     firstVertId++;
  }
}

void
TriFacet::PrintTriFacetTopo(FILE* file=stderr, int& TriFacetId, int& firstVertId)
{
  fprintf(file," %6d  %6e  %6e  %6e\n",TriFacetId,firstVertId,firstVertId+1,firstVertId+2);
  TriFacetId++;
  firstVertId += 3;
}
*/

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
TriFacet::ComputeApproxArea(CoordSet &cs)
{
  double M1[3], M2[3], M3[3];
  double *m;
  m = &LocalCoordOnFaceEl[0][0];
  (*FaceEl).LocalToGlobalCoord(M1, m, cs);
  m = &LocalCoordOnFaceEl[1][0];
  (*FaceEl).LocalToGlobalCoord(M2, m, cs);
  m = &LocalCoordOnFaceEl[2][0];
  (*FaceEl).LocalToGlobalCoord(M3, m, cs);

  Vector V1(M1,3), V2(M2,3), V3(M3,3);
  Vector V12 = V2 - V1;  
  Vector V13 = V3 - V1; 
  Vector W = V12.cross(V13);

  double S = 0.5*W.magnitude();

  SetArea(S);
  return S;
}*/


// -----------------------------------------------------------------------------------------------------
//                         REF. TRIANGULAR FACET -> REF. FACE El. MAPPING METHODS 
// -----------------------------------------------------------------------------------------------------
/*void
TriFacet::MappingShapeFctAndDerivative(double* Shape, double* dShape, double* m)
{
  if(Shape==0)  Shape  = new double[3];
  if(dShape==0) dShape = new double[2][3];

  Shape[0] = m[0]; Shape[1] = m[1]; Shape[2] = 1.-m[0]-m[1];
  dShape[0][0] =  1.; dShape[0][1] =  0.; 
  dShape[1][0] =  0.; dShape[1][1] =  1.; 
  dShape[2][0] = -1.; dShape[2][1] = -1.; 
}*/

// Return jacobian of the mapping ref. TriFacet -> ref. face el.
template<class Scalar>
Scalar
TriFacet<Scalar>::MappingJacobian()
{
  // J = 2.Area = ||12 x 13|| 
  Scalar X[3] = {LocalCoordOnFaceEl[0][0],LocalCoordOnFaceEl[1][0],LocalCoordOnFaceEl[2][0]};
  Scalar Y[3] = {LocalCoordOnFaceEl[0][1],LocalCoordOnFaceEl[1][1],LocalCoordOnFaceEl[2][1]};

  //cerr << " SENS A VERIFIER !!!" << endl;
  Scalar J = (X[1]-X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]);

  using std::abs;
  //if(J < 0) std::cerr << "here in TriFacet<Scalar>::MappingJacobian, J < 0\n";
  return abs(J);// to avoid swapping the vertex to have a positive jacobian
                // this is OK because we just use it as the differential area
                // in the integration formula
}

// Map the given point m in the ref. Trifacet to the associated point in the ref. face el. 
template<class Scalar>
void
TriFacet<Scalar>::LocalToLocalCoordOnFaceEl(const double* m, Scalar* mOnFaceEl)
{
  double r = m[0];
  double s = m[1];
  double t = 1 - r - s;

  mOnFaceEl[0] = t*LocalCoordOnFaceEl[0][0] + r*LocalCoordOnFaceEl[1][0] + s*LocalCoordOnFaceEl[2][0];
  mOnFaceEl[1] = t*LocalCoordOnFaceEl[0][1] + r*LocalCoordOnFaceEl[1][1] + s*LocalCoordOnFaceEl[2][1];
}

// Return the jacobian on the real face el. at the point associated to 
// the given point m in the ref. TriFacet 
template<class Scalar>
template<typename CoordSetT>
Scalar
TriFacet<Scalar>::GetJacobianOnFaceEl(const Scalar* m, CoordSetT &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping TriFacet -> ref. face el.)
  Scalar JTriFacetMapping = MappingJacobian();
   
  Scalar mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  Scalar JFaceElMapping = FaceEl->GetJacobian(mOnFaceEl, cs);

  return JFaceElMapping*JTriFacetMapping;
}

// Return the isoparametric mapping normal and jacobian on the real face el. at the point 
// associated to the given point m in the ref. TriFacet
template<class Scalar>
template<typename CoordSetT>
Scalar
TriFacet<Scalar>::GetIsoParamMappingNormalAndJacobianOnFaceEl(Scalar* Normal, const Scalar* m, CoordSetT &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping TriFacet -> ref. face el.)
  Scalar JTriFacetMapping = MappingJacobian();

  Scalar mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  Scalar JFaceElMapping = FaceEl->GetIsoParamMappingNormalAndJacobian(Normal, mOnFaceEl, cs);

  return JFaceElMapping*JTriFacetMapping;
}

// ------------------------------------------------------------------------------------------------------
//                         INTEGRATION OF SHAPE FUNCTIONS PRODUCT METHODS 
// ------------------------------------------------------------------------------------------------------
template<>
inline FullM
TriFacet<double>::IntegrateShapeFctProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSet &cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// NOTE: 
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet 
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFaceElem = FaceElem->nNodes();

   //cerr << "In TriFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFaceElem = " << nShapeFctOnFaceElem << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,nShapeFctOnFaceElem);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double)); 
   double* ShapeOnFaceElem  = (double*) dbg_alloca(nShapeFctOnFaceElem*sizeof(double)); 
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;

   int igp, i, j;
   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;
   
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetJacobianOnFaceEl(m, cs); 
     //cerr << " # dA = " << dA << endl;
 
     // Get local coord. on each face element 
     // -> for the mortar elem. (see the NOTE section) 
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;
 
     // -> for the element supporting the given triangular facet 
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl); 
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;
     
     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++)
     //  cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;
 
     // -> for the element supporting the given triangular facet
     FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++)
     //  cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;
      
     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFaceElem;j++)
         MatShapeFctProd[i][j] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j];
   }  
 
  //MatShapeFctProd.print("M[TriFacet]="); 
  return MatShapeFctProd; 
}

#if (MAX_MORTAR_DERIVATIVES > 0)
template<>
inline FullM
TriFacet<ActiveDouble>::IntegrateShapeFctProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSet &cs, int ngp)
{
  std::cerr << "TriFacet<ActiveDouble>::IntegrateShapeFctProduct is not implemented\n";
  return FullM(0);
}
#endif

template<>
inline FullM
TriFacet<double>::IntegrateNormalShapeFctProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSet& cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current TriFacet][Mortar(i).FriendFacet.Shape(j).normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO 
//                                       THE GIVEN (FRIEND) TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFaceElem = FaceElem->nNodes();
   //cerr << "In TriFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFaceElem = " << nShapeFctOnFaceElem << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nShapeFctOnFaceElem);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double));
   double* ShapeOnFaceElem  = (double*) dbg_alloca(nShapeFctOnFaceElem*sizeof(double));
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;
   double Normal[3];
   int igp, i, j;

   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;
  
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetIsoParamMappingNormalAndJacobianOnFaceEl(Normal, m, cs);
     //cerr << " # dA = " << dA << endl;
     //cerr << " # normal = " << Normal[0] <<" "<< Normal[1] <<" "<< Normal[2] << endl;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++)
     //  cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;

     // -> for the element supporting the given triangular facet
     FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++)
     //  cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;

     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFaceElem;j++){
         MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[0];
         MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[1];
         MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[2];
       } 
  }

  //MatShapeFctProd.print("M[TriFacet]=");
  return MatShapeFctProd;
}

#if (MAX_MORTAR_DERIVATIVES > 0)
template<>
inline FullM
TriFacet<ActiveDouble>::IntegrateNormalShapeFctProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSet& cs, int ngp)
{
  std::cerr << "TriFacet<ActiveDouble>::IntegrateNormalShapeFctProduct is not implemented\n";
  return FullM(0);
}
#endif

template<>
inline FullM
TriFacet<double>::IntegrateGradNormalShapeFctProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSet& cs,
                                                     CoordSet& cs1, TriFacet& FriendFacet2, CoordSet& cs2, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the gradient of the shape functions defined by
// the given MortarElement and the shape functions of the element associated to the given (Friend) triangular
// facet times the normal (see Notes (2)).
// -> Mij = Intg[current TriFacet][Mortar(i).FriendFacet.Shape(j).normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO 
//                                       THE GIVEN (FRIEND) TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
//     (3) cs1 is the CoordSet of the FriendFacet, cs2 is the CoordSet of FriendFacet2
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FriendFaceEl = FriendFacet.GetPtrFaceEl();
   // Get ptr to the face element supporting the second given triangular facet
   FaceElement* FriendFaceEl2 = FriendFacet2.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
   //cerr << "In TriFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFriendFaceEl = " << nShapeFctOnFriendFaceEl << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct, 3*nShapeFctOnFriendFaceEl);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double));
   double* ShapeOnFriendFaceEl  = (double*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(double));
   double m[2], mOnMortarEl[2], mOnFriendFaceEl[2];
   double r, s, t, weight;
   double JNormal[3];

   int nShapeFctOnFriendFaceEl2 = FriendFaceEl->nNodes();
   double* ShapeOnFriendFaceEl2  = (double*) dbg_alloca(nShapeFctOnFriendFaceEl2*sizeof(double));
   double mOnFriendFaceEl2[2], MOnFriendFaceEl[3], MOnFriendFaceEl2[3], dJNormal[12][3];

   int igp, i, j;
  
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);
     //cerr << " # mOnFriendFaceEl:   x = " << mOnFriendFaceEl[0] << ", y = " << mOnFriendFaceEl[1] << endl;

     // -> for the element supporting the second given triangular facet
     FriendFacet2.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl2);
     //cerr << " # mOnFriendFaceEl:   x = " << mOnFriendFaceEl[0] << ", y = " << mOnFriendFaceEl[1] << endl;

     // Jacobian and JNormal on the face element supporting the CURRENT triangular facet
     double dA = MappingJacobian(); 
     FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
     //cerr << " # dA = " << dA << endl;
     //cerr << " # J*normal = " << JNormal[0] <<" "<< JNormal[1] <<" "<< JNormal[2] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++) cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;

     // -> for the element supporting the given triangular facet
     FriendFaceEl->GetShapeFctVal(ShapeOnFriendFaceEl, mOnFriendFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++) cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;

     // -> for the element supporting the second given triangular facet
     FriendFaceEl2->GetShapeFctVal(ShapeOnFriendFaceEl2, mOnFriendFaceEl2);
     //for(j=0;j<nShapeFctOnFriendFaceEl2;j++) cerr << " ShapeOnFriendFaceEl2[" << j << "] = " << ShapeOnFaceElem2[j] << endl;

     // Get the derivative of J*Normal
     FaceEl->GetdJNormal(dJNormal, mOnMortarEl, cs);
     //for(int i=0; i<12; ++i) cerr << " # dJNormal[" << i << "] = " << dJNormal[i][0] <<" "<< dJNormal[i][1] <<" "<< dJNormal[i][2] << endl;

     // Get global coord. on the element supporting the given triangular facet
     FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs1);
     //cerr << " # MOnFriendFaceEl:   x = " << MOnFriendFaceEl[0] << ", y = " << MOnFriendFaceEl[1] << ", z = " << MOnFriendFaceEl[2] << endl;

     // Get global coord. on the element supporting the given other triangular facet
     FriendFaceEl2->LocalToGlobalCoord(MOnFriendFaceEl2, mOnFriendFaceEl2, cs2);
     //cerr << " # MOnFriendFaceEl:   x = " << MOnFriendFaceEl[0] << ", y = " << MOnFriendFaceEl[1] << ", z = " << MOnFriendFaceEl[2] << endl;

     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFriendFaceEl;j++){
         MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*ShapeOnFriendFaceEl[j]*JNormal[0];
         MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*ShapeOnFriendFaceEl[j]*JNormal[1];
         MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*ShapeOnFriendFaceEl[j]*JNormal[2];
         //cerr << "i = " << i << ", j = " << j << ", MatShapeFctProd[i] #1 = " << MatShapeFctProd[i][3*j  ] << " " << MatShapeFctProd[i][3*j+1]
         //     << " " << MatShapeFctProd[i][3*j+2] << endl;
         MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*(MOnFriendFaceEl[0]*dJNormal[3*j  ][0]
                                      +MOnFriendFaceEl[1]*dJNormal[3*j  ][1]+MOnFriendFaceEl[2]*dJNormal[3*j  ][2]);
         MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*(MOnFriendFaceEl[0]*dJNormal[3*j+1][0]
                                      +MOnFriendFaceEl[1]*dJNormal[3*j+1][1]+MOnFriendFaceEl[2]*dJNormal[3*j+1][2]);
         MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*(MOnFriendFaceEl[0]*dJNormal[3*j+2][0]
                                      +MOnFriendFaceEl[1]*dJNormal[3*j+2][1]+MOnFriendFaceEl[2]*dJNormal[3*j+2][2]);
         MatShapeFctProd[i][3*j  ] -= weight*dA*MortarShape[i]*(MOnFriendFaceEl2[0]*dJNormal[3*j  ][0]
                                      +MOnFriendFaceEl2[1]*dJNormal[3*j  ][1]+MOnFriendFaceEl2[2]*dJNormal[3*j  ][2]);
         MatShapeFctProd[i][3*j+1] -= weight*dA*MortarShape[i]*(MOnFriendFaceEl2[0]*dJNormal[3*j+1][0]
                                      +MOnFriendFaceEl2[1]*dJNormal[3*j+1][1]+MOnFriendFaceEl2[2]*dJNormal[3*j+1][2]);
         MatShapeFctProd[i][3*j+2] -= weight*dA*MortarShape[i]*(MOnFriendFaceEl2[0]*dJNormal[3*j+2][0]
                                      +MOnFriendFaceEl2[1]*dJNormal[3*j+2][1]+MOnFriendFaceEl2[2]*dJNormal[3*j+2][2]);
         //cerr << "i = " << i << ", j = " << j << ", MatShapeFctProd[i] #2 = " << MatShapeFctProd[i][3*j  ] << " " << MatShapeFctProd[i][3*j+1]
         //     << " " << MatShapeFctProd[i][3*j+2] << endl;
       } 
  }

  //MatShapeFctProd.print("M[TriFacet]=");
  return MatShapeFctProd;
}

/* old version
// Compute normal "geometrical" gap 
template<class Scalar>
template<typename CoordSetT>
GenVector<Scalar>
TriFacet<Scalar>::IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, TriFacet& FriendFacet, CoordSetT& cs, CoordSetT& cs1, int ngp,
                                                double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current TriFacet][Mortar(i).FriendFacet.Shape(j).normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO
//                                       THE GIVEN (FRIEND) TRIFACET)
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FriendFaceEl = FriendFacet.GetPtrFaceEl();
   int nMortarShapeFct   = MortarEl->nNodes();

   GenVector<Scalar> NormalGeoGaps(nMortarShapeFct,Scalar(0.0));

   Scalar* MortarShape = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
   Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
   double m[2], r, s, t, weight;
   Scalar J, JNormal[3], MOnFriendFaceEl[3];
   int igp, i;

   Scalar dA = MappingJacobian();
   if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
     FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
     if(offset != 0) J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
   }

   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl); 
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl); 
     //cerr << " # mOnFriendFaceEl:   x = " << mOnFriendFaceEl[0] << ", y = " << mOnFriendFaceEl[1] << endl;

     // Jacobian and J*Normal on the face element supporting the CURRENT triangular facet
     if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
       FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
       if(offset != 0) J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
     }
     //cerr << " # dA = " << dA << ", J*Normal = " << JNormal[0] <<" "<< JNormal[1] <<" "<< JNormal[2] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //cerr << " # MortarShape = "; for(i=0; i<nMortarShapeFct; ++i) cerr << MortarShape[i] << " "; cerr << endl;

     // Get global coord. on the element supporting the given triangular facet
     FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs1); // PJSA changed cs to cs1 (Friend's CoordSet)
     //cerr << " # MOnFriendFaceEl:   x = " << MOnFriendFaceEl[0] << ", y = " << MOnFriendFaceEl[1] << ", z = " << MOnFriendFaceEl[2] << endl;

     // Shape fcts product & integration
     Scalar sdot = weight*dA*(MOnFriendFaceEl[0]*JNormal[0]+MOnFriendFaceEl[1]*JNormal[1]+MOnFriendFaceEl[2]*JNormal[2]);
     //cerr << " # sdot = " << sdot << endl;
     for(i=0;i<nMortarShapeFct;i++) {
       NormalGeoGaps[i] += sdot*MortarShape[i];
       if(offset != 0) NormalGeoGaps[i] += weight*dA*J*MortarShape[i]*offset;
     }
  }

  //NormalGeoGaps.print("  NormalGeoGaps[TriFacet]=");
  return NormalGeoGaps;
}
*/

// Compute normal "geometrical" gap  in one pass by integrating over A and B simultaneously
template<class Scalar>
template<typename CoordSetT>
GenVector<Scalar>
TriFacet<Scalar>::IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, TriFacet& FriendFacetB, TriFacet& FriendFacetC,
                                                CoordSetT& cs, CoordSetT& csB, CoordSetT& csC, int ngp,
                                                double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current TriFacet][Mortar(i).FriendFacet.Shape(j).normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO
//                                       THE GIVEN (FRIEND) TRIFACET)
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facets
   FaceElement* FriendFaceElB = FriendFacetB.GetPtrFaceEl();
   FaceElement* FriendFaceElC = FriendFacetC.GetPtrFaceEl();
   int nMortarShapeFct   = MortarEl->nNodes();

   GenVector<Scalar> NormalGeoGaps(nMortarShapeFct,Scalar(0.0));

   Scalar* MortarShape = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
   Scalar mOnMortarEl[2], mOnFriendFaceElB[2], mOnFriendFaceElC[2];
   double m[2], r, s, t, weight;
   Scalar J, JNormal[3], MOnFriendFaceElB[3], MOnFriendFaceElC[3];
   int igp, i;

   Scalar dA = MappingJacobian();
   if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
     FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
     if(offset != 0) J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
   }

   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl); 
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facets
     FriendFacetB.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceElB); 
     FriendFacetC.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceElC);
     //cerr << " # mOnFriendFaceElB:   x = " << mOnFriendFaceElB[0] << ", y = " << mOnFriendFaceElB[1] << endl;
     //cerr << " # mOnFriendFaceElC:   x = " << mOnFriendFaceElC[0] << ", y = " << mOnFriendFaceElC[1] << endl;

     // Jacobian and J*Normal on the face element supporting the CURRENT triangular facet
     if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
       FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
       if(offset != 0) J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
     }
     //cerr << " # dA = " << dA << ", J*Normal = " << JNormal[0] <<" "<< JNormal[1] <<" "<< JNormal[2] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //cerr << " # MortarShape = "; for(i=0; i<nMortarShapeFct; ++i) cerr << MortarShape[i] << " "; cerr << endl;

     // Get global coord. on the element supporting the given triangular facets
     FriendFaceElB->LocalToGlobalCoord(MOnFriendFaceElB, mOnFriendFaceElB, csB);
     FriendFaceElC->LocalToGlobalCoord(MOnFriendFaceElC, mOnFriendFaceElC, csC);
     //cerr << " # MOnFriendFaceElB:   x = " << MOnFriendFaceElB[0] << ", y = " << MOnFriendFaceElB[1] << ", z = " << MOnFriendFaceElB[2] << endl;
     //cerr << " # MOnFriendFaceElC:   x = " << MOnFriendFaceElC[0] << ", y = " << MOnFriendFaceElC[1] << ", z = " << MOnFriendFaceElC[2] << endl;

     // Shape fcts product & integration
     Scalar sdot = weight*dA*( (MOnFriendFaceElB[0]-MOnFriendFaceElC[0])*JNormal[0]
                              +(MOnFriendFaceElB[1]-MOnFriendFaceElC[1])*JNormal[1]
                              +(MOnFriendFaceElB[2]-MOnFriendFaceElC[2])*JNormal[2] );
     //cerr << " # sdot = " << sdot << endl;
     for(i=0;i<nMortarShapeFct;i++) {
       NormalGeoGaps[i] += sdot*MortarShape[i];
       if(offset != 0) NormalGeoGaps[i] += weight*dA*J*MortarShape[i]*offset;
     }
  }

  //NormalGeoGaps.print("  NormalGeoGaps[TriFacet]=");
  return NormalGeoGaps;
}
