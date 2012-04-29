// Std C/C++ headers
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
//#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/FFIPolygon.d/TriFacet.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
FFIPolygon<Scalar>::FFIPolygon()
: Facets()
{
  Initialize();
}

template<class Scalar>
FFIPolygon<Scalar>::FFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data)
: Facets()
{
  Initialize();
  SetFFIPolygon(MasterFaceEl, SlaveFaceEl, nVert, ACME_FFI_Data);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
FFIPolygon<Scalar>::~FFIPolygon()
{
#ifdef HB_ACME_FFI_DEBUG
  if(VertexLlCoordOnSFaceEl) { delete VertexLlCoordOnSFaceEl; }
  if(VertexLlCoordOnMFaceEl) { delete VertexLlCoordOnMFaceEl; }
#endif
  if(dM) delete [] dM;
  if(dN) delete [] dN;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void 
FFIPolygon<Scalar>::Initialize()
{
  Area      = 0.0;
  nVertices = 0;
  MasterFace= 0;
  SlaveFace = 0;

#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = 0;
  VertexLlCoordOnMFaceEl = 0;
#endif

  dM = 0;
  dN = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
FFIPolygon<Scalar>::SetPtrMasterFace(FaceElement* PtrMasterFace) { MasterFace = PtrMasterFace; }

template<class Scalar>
void
FFIPolygon<Scalar>::SetPtrSlaveFace(FaceElement* PtrSlaveFace) { SlaveFace = PtrSlaveFace; }



template<class Scalar>
void
FFIPolygon<Scalar>::SetFFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, 
                                  int nVert, double* ACME_FFI_Data)
{
  nVertices  = nVert;
  MasterFace = MasterFaceEl;
  SlaveFace  = SlaveFaceEl;

  // Create & fill polygon triangularization
  int LocalCoordOffset = 2*nVertices+1;
 
  double* ACME_FFI_LocalCoordData = &ACME_FFI_Data[LocalCoordOffset];
  Scalar* ACME_FFI_LocalCoordData_copy = new Scalar[nVertices*4];

  for(int i=0; i<nVertices; ++i) {
    ACME_FFI_LocalCoordData_copy[4*i  ] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+1] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+2] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+3] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData += 4*MAX_FFI_DERIVATIVES;
  }
  CreateTriangularization(ACME_FFI_LocalCoordData_copy);
  delete [] ACME_FFI_LocalCoordData_copy;
}

#if (MAX_MORTAR_DERIVATIVES > 0)
template<>
inline void
FFIPolygon<MadDouble>::SetFFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl,
                                     int nVert, double* ACME_FFI_Data)
{
  nVertices  = nVert;
  MasterFace = MasterFaceEl;
  SlaveFace  = SlaveFaceEl;

  // Create & fill polygon triangularization
  int LocalCoordOffset = 2*nVertices+1;
  ACME_FFI_LocalCoordData = &ACME_FFI_Data[LocalCoordOffset];
}

template<>
inline void
FFIPolygon<MadDouble>::PrepLocalCoordData()
{
  MadDouble* ACME_FFI_LocalCoordData_copy = new MadDouble[nVertices*4];

  for(int i=0; i<nVertices; ++i) {
    ACME_FFI_LocalCoordData_copy[4*i  ] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+1] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+2] = *ACME_FFI_LocalCoordData++;
    ACME_FFI_LocalCoordData_copy[4*i+3] = *ACME_FFI_LocalCoordData++;
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+0].fastAccessDx(j) = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+1].fastAccessDx(j) = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+2].fastAccessDx(j) = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+3].fastAccessDx(j) = *ACME_FFI_LocalCoordData++;
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+0].derivatives()[j] = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+1].derivatives()[j] = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+2].derivatives()[j] = *ACME_FFI_LocalCoordData++;
    for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) ACME_FFI_LocalCoordData_copy[4*i+3].derivatives()[j] = *ACME_FFI_LocalCoordData++;
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
    for(int k=0; k<4; ++k)
      for(int j=0; j<MAX_FFI_DERIVATIVES; ++j)
       ACME_FFI_LocalCoordData++; // TODO need to use the first derivatives of the FFI, but how?
    for(int k=0; k<4; ++k) {
      for(int j=0; j<MAX_FFI_DERIVATIVES; ++j)
        for(int l=0; l<=j; ++l) {
          ACME_FFI_LocalCoordData++; // TODO need to use the second derivatives of the FFI, but how?
        }
    }
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
    for(int k=0; k<4; ++k)
      for(int j=0; j<MAX_FFI_DERIVATIVES; ++j) {
        ACME_FFI_LocalCoordData_copy[4*i+k].derivatives()[j].value() = *ACME_FFI_LocalCoordData;
        ACME_FFI_LocalCoordData_copy[4*i+k].value().derivatives()[j] = *ACME_FFI_LocalCoordData++; // NOTE!!!
      }
    for(int k=0; k<4; ++k) {
      for(int j=0; j<MAX_FFI_DERIVATIVES; ++j)
        for(int l=0; l<=j; ++l) {
          ACME_FFI_LocalCoordData_copy[4*i+k].derivatives()[j].derivatives()[l] = *ACME_FFI_LocalCoordData;
          ACME_FFI_LocalCoordData_copy[4*i+k].derivatives()[l].derivatives()[j] = *ACME_FFI_LocalCoordData++;
        }
    }
#else
    // TODO: MORTAR_AUTO_DIFF_SACADO_RAD_FAD
    ACME_FFI_LocalCoordData += (4*MAX_FFI_DERIVATIVES+4*MAX_FFI_SECOND_DERIVATIVES);
#endif
  }

  CreateTriangularization(ACME_FFI_LocalCoordData_copy);
  delete [] ACME_FFI_LocalCoordData_copy;
}
#endif

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
FFIPolygon<Scalar>::Print()
{
  fprintf(stderr," * -------------------------------- \n");
  fprintf(stderr," * FFIPolygon data:\n");
  fprintf(stderr," + number of vertices: %d\n", GetnVertices()); 
  fprintf(stderr," + number of TriFacet: %d\n", GetnFacets()); 
  fprintf(stderr," + ptr to master elem: %p\n", MasterFace);
  if(MasterFace) MasterFace->print();
  fprintf(stderr," + ptr to slave elem : %p\n", SlaveFace);
  if(SlaveFace) SlaveFace->print();
  fprintf(stderr," + TriFacets data on master face:\n");  
  for(int iTri=0; iTri<GetnFacets(); iTri++){
    fprintf(stderr," -> TriFacet %d (on master) \n",iTri+1);
    MasterFacet(Facets,iTri).Print();    
  } 
  fprintf(stderr," + TriFacets data on slave face:\n");  
  for(int iTri=0; iTri<GetnFacets(); iTri++){
    fprintf(stderr," -> TriFacet %d (on slave) \n",iTri+1);
    SlaveFacet(Facets,iTri).Print();    
  } 
  fprintf(stderr," * -------------------------------- \n");
}

template<class Scalar>
void
FFIPolygon<Scalar>::PrintM()
{
  GetPtrSlaveFace()->print(); 
  GetPtrSlaveFace()->print(); 
  M.print("M[FFI] = ");
}

template<class Scalar>
void
FFIPolygon<Scalar>::PrintN()
{
  GetPtrSlaveFace()->print(); 
  GetPtrMasterFace()->print(); 
  N.print("N[FFI] = ");
}

#ifdef HB_ACME_FFI_DEBUG
template<class Scalar>
void 
FFIPolygon<Scalar>::PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnSFaceEl[iVert][0],VertexLlCoordOnSFaceEl[iVert][1]};
    SlaveFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

template<class Scalar>
void
FFIPolygon<Scalar>::PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnMFaceEl[iVert][0],VertexLlCoordOnMFaceEl[iVert][1]};
    MasterFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

template<class Scalar>
void
FFIPolygon<Scalar>::PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode)
{
   int firstVertId = VertOffset;
   switch(nVertices) {
     case 4 : { 
       fprintf(file," %6d  %3d  %6d  %6d  %6d  %6d\n",EdgeOffset,188,VertOffset,VertOffset+1,VertOffset+2,VertOffset+3);
       VertOffset+=4; EdgeOffset++;
     } break;
     case 3 : {
       fprintf(file," %6d  %3d  %6d  %6d  %6d\n",EdgeOffset,108,VertOffset,VertOffset+1,VertOffset+2);
       VertOffset+=3; EdgeOffset++;
     } break;
     default : {
       for(int iVert=0; iVert<nVertices-1; iVert++){
         fprintf(file," %6d  %3d  %6d  %6d\n",EdgeOffset,elCode,VertOffset,VertOffset+1);
         VertOffset++; EdgeOffset++;
       }
       fprintf(file," %6d  %3d  %6d  %6d\n",EdgeOffset,elCode,VertOffset,firstVertId);
       VertOffset++; EdgeOffset++;
     } break;
   }
}
#endif

// -----------------------------------------------------------------------------------------------------
//                                       MEMORY ALLOCATION/DESALOCATION  METHODS 
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
//                                       TRIANGULARIZATION METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
void
FFIPolygon<Scalar>::CreateTriangularization(Scalar* ACME_FFI_LocalCoordData)
{
  // Allocate space for Facets
  // By default assume as many triangular facet as polygon vertices
  // -> can be overwritten in for optimization
#ifdef MORTAR_AUTO_DIFF_SACADO_RAD_FAD // workaround unresolved issue with Sacado::Rad and stl::vector 
  Facet_pair_t dummy;
  for(int i=0; i<nVertices; ++i) Facets.push_back(dummy);
#else
  Facets.assign(nVertices, Facet_pair_t(TriFacet<Scalar>(), TriFacet<Scalar>()));
#endif

  // Create master & slave polygon centroids
  Scalar MasterCentroid[2] = {0.0, 0.0};
  Scalar SlaveCentroid[2]  = {0.0, 0.0};
#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = new Scalar[nVertices][2];
  VertexLlCoordOnMFaceEl = new Scalar[nVertices][2];
#endif
  for(int iVert=0; iVert<nVertices; ++iVert) {
     int offset = 4*iVert;
#ifdef HB_ACME_FFI_DEBUG
     VertexLlCoordOnMFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset  ];
     VertexLlCoordOnMFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+1];
     VertexLlCoordOnSFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset+2];
     VertexLlCoordOnSFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+3];
#endif 
     MasterCentroid[0] += ACME_FFI_LocalCoordData[offset  ]; 
     MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+1]; 
     SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset+2]; 
     SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+3]; 
  }
  SlaveCentroid[0]  /= Scalar(nVertices); SlaveCentroid[1]  /= Scalar(nVertices);
  MasterCentroid[0] /= Scalar(nVertices); MasterCentroid[1] /= Scalar(nVertices);

  //std::cerr << " SlaveCentroid : x = " << SlaveCentroid[0] << ", y = " << SlaveCentroid[1] << std::endl;
  //std::cerr << " MasterCentroid : x = " << MasterCentroid[0] << ", y = " << MasterCentroid[1] << std::endl;

  // Create master & slave triangularization
  for(size_t i=0, ii=GetnFacets(), offsetm1=0; i<ii; i++, offsetm1+=4) {
     size_t offsetm2 = (i==ii-1) ? 0 : offsetm1+4;
     //HB: note that ACME Master <=> our Slave, ... 
     MasterFacet(Facets, i).SetTriFacet(MasterFace, MasterCentroid, &ACME_FFI_LocalCoordData[offsetm1],
                                                                    &ACME_FFI_LocalCoordData[offsetm2]);
     SlaveFacet(Facets , i).SetTriFacet(SlaveFace,  SlaveCentroid,  &ACME_FFI_LocalCoordData[offsetm1+2],
                                                                    &ACME_FFI_LocalCoordData[offsetm2+2]);
  } 
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FFIPolygon<Scalar>::ComputeArea()
{
   double Area = 0.0;
   for(size_t i=0, ii=GetnFacets(); i<ii; i++) {
      Area += TriFacet[i]->ComputeArea();
   }
   return Area;
}*/

// -----------------------------------------------------------------------------------------------------
//                         INTEGRATION OF SHAPE FUNCTIONS PRODUCT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnMaster_MasterShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the MASTER SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

   FullM MatShapeFctProd(nMortarShapeFct,nMasterShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += MasterFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                        cs, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnMaster_SlaveShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the MASTER SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nSlaveShapeFct  = SlaveFace->nNodes();

   FullM MatShapeFctProd(nMortarShapeFct,nSlaveShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += MasterFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                        cs, ngp);
   } 

   return MatShapeFctProd;
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

   //cerr << "In FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nMasterShapeFct = " << nMasterShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;
   
   FullM MatShapeFctProd(nMortarShapeFct,nMasterShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                       cs, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_SlaveShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nSlaveShapeFct  = SlaveFace->nNodes();
   
   //cerr << "In FFIPolygon<Scalar>::IntegrateOnSlave_SlaveShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nSlaveShapeFct  = " << nSlaveShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,nSlaveShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                       cs, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
void 
FFIPolygon<Scalar>::ComputeM(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Slave shape fcts product mass-like matrix M
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_SlaveShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_SlaveShapeFctProduct for more details
// *****************************************************************************************************
{
  M = IntegrateOnSlave_SlaveShapeFctProduct(MortarEl, cs, ngp);
  //M.print();
}

template<class Scalar>
void 
FFIPolygon<Scalar>::ComputeN(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterShapeFctProduct(MortarEl, cs, ngp);
  //N.print();
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_MasterNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the normal and the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element 
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

   //cerr << "In FFIPolygon::IntegrateOnSlave_MasterNormalShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nMasterShapeFct = " << nMasterShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nMasterShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateNormalShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                             cs, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the normal and the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nSlaveShapeFct  = SlaveFace->nNodes();                                                                                                                                          
   //cerr << "In FFIPolygon::IntegrateOnSlave_SlaveNormalShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nMasterShapeFct = " << nMasterShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nSlaveShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateNormalShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                             cs, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
void
FFIPolygon<Scalar>::ComputeNormalM(MortarElement* MortarEl, CoordSet &cs, int ngp)
{
  M = IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarEl, cs, ngp);
}

template<class Scalar>
void
FFIPolygon<Scalar>::ComputeNormalN(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygon<Scalar>::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterNormalShapeFctProduct(MortarEl, cs, ngp);
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_MasterGradNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the contribution to the gradient of the gap function term due to the
// product of the normal and the shape functions defined by the given MortarElement and the shape
// functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

   //cerr << "In FFIPolygon::IntegrateOnSlave_MasterGradNormalShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nMasterShapeFct = " << nMasterShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nMasterShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateGradNormalShapeFctProduct(MortarEl, MasterFacet(Facets, i), SlaveCoords,
                                                                                 MasterCoords, MasterFacet(Facets, i), MasterCoords, ngp);
   }

   return MatShapeFctProd;
}

template<class Scalar>
FullM
FFIPolygon<Scalar>::IntegrateOnSlave_SlaveGradNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the contribution to the gradient of the gap function term due to the
// product of the normal and the shape functions defined by the given MortarElement and the shape
// functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nSlaveShapeFct  = SlaveFace->nNodes();

   //cerr << "In FFIPolygon::IntegrateOnSlave_SlaveGradNormalShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct = " << nMortarShapeFct << endl;
   //cerr << " -> nMasterShapeFct = " << nMasterShapeFct << endl;
   //cerr << " -> ngp             = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nSlaveShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += SlaveFacet(Facets, i).IntegrateGradNormalShapeFctProduct(MortarEl, SlaveFacet(Facets, i), SlaveCoords,
                                                                                 SlaveCoords, MasterFacet(Facets, i), MasterCoords, ngp);
   }

   return MatShapeFctProd;
}

template<>
inline void
FFIPolygon<double>::ComputeGradNormalM(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
  M = IntegrateOnSlave_SlaveGradNormalShapeFctProduct(MortarEl, SlaveCoords, MasterCoords, ngp);
}

template<>
inline void
FFIPolygon<double>::ComputeGradNormalN(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
  N = IntegrateOnSlave_MasterGradNormalShapeFctProduct(MortarEl, SlaveCoords, MasterCoords, ngp);
}

template<>
inline void
FFIPolygon<double>::ComputeNormalGeoGap(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp,
                                        double offset)
{
   int nMortarShapeFct = MortarEl->nNodes();

   NormalGeoGaps.reset(nMortarShapeFct,0.0);

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
/*   // old version
     NormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, SlaveFacet(Facets, i), SlaveCoords,
                                                                          SlaveCoords, ngp, offset);
     NormalGeoGaps -= SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, MasterFacet(Facets, i), SlaveCoords,
                                                                          MasterCoords, ngp);*/
     // new version
     NormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, SlaveFacet(Facets, i), MasterFacet(Facets, i),
                                                                          SlaveCoords, SlaveCoords, MasterCoords, ngp, offset);
   }
   //NormalGeoGaps.print(" NormalGeoGaps[FFIPolygon]");
}

#if (MAX_MORTAR_DERIVATIVES > 0)
template<>
inline void
FFIPolygon<MadDouble>::ComputeGradNormalM(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
   // M is computed in FFIPolygon<MadDouble>::ComputeNormalGeoGap
}

template<>
inline void
FFIPolygon<MadDouble>::ComputeGradNormalN(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
   // N is computed in FFIPolygon<MadDouble>::ComputeNormalGeoGap
}

template<>
inline void
FFIPolygon<MadDouble>::ComputeNormalGeoGap(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp,
                                           double offset)
{
   PrepLocalCoordData();
   int nMortarShapeFct = MortarEl->nNodes();

   GenVector<MadDouble> MadNormalGeoGaps(nMortarShapeFct,MadDouble(0.0));

   MadCoordSet MadMasterCoords, MadSlaveCoords;
   for (int i=0; i<MasterFace->nNodes(); ++i) {
     Node& node = MasterCoords.getNode(MasterFace->GetNode(i));
     MadNode *mynode = new MadNode(node.x,node.y,node.z);
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
     mynode->x.diff(3*i+0, MAX_MORTAR_DERIVATIVES);
     mynode->y.diff(3*i+1, MAX_MORTAR_DERIVATIVES);
     mynode->z.diff(3*i+2, MAX_MORTAR_DERIVATIVES);
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
     mynode->x = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, 3*i+0, node.x);
     mynode->y = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, 3*i+1, node.y);
     mynode->z = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, 3*i+2, node.z);
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
     mynode->x.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, 3*i+0);
     mynode->y.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, 3*i+1);
     mynode->z.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, 3*i+2);
#else // MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
     mynode->x = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.x,MAX_MORTAR_DERIVATIVES,3*i+0), MAX_MORTAR_DERIVATIVES, 3*i+0);
     mynode->y = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.y,MAX_MORTAR_DERIVATIVES,3*i+1), MAX_MORTAR_DERIVATIVES, 3*i+1);
     mynode->z = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.z,MAX_MORTAR_DERIVATIVES,3*i+2), MAX_MORTAR_DERIVATIVES, 3*i+2);
#endif
     MadMasterCoords[MasterFace->GetNode(i)] = mynode;
   }
   for (int i=0; i<SlaveFace->nNodes(); ++i) {
     Node& node = SlaveCoords.getNode(SlaveFace->GetNode(i));
     MadNode *mynode = new MadNode(node.x,node.y,node.z);
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
     mynode->x.diff(MAX_MORTAR_DERIVATIVES/2+3*i+0, MAX_MORTAR_DERIVATIVES);
     mynode->y.diff(MAX_MORTAR_DERIVATIVES/2+3*i+1, MAX_MORTAR_DERIVATIVES);
     mynode->z.diff(MAX_MORTAR_DERIVATIVES/2+3*i+2, MAX_MORTAR_DERIVATIVES);
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
     mynode->x = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+0, node.x);
     mynode->y = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+1, node.y);
     mynode->z = Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES>(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+2, node.z);
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
     mynode->x.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+0);
     mynode->y.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+1);
     mynode->z.derivatives() = DerivativeType::Unit(MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+2);
#else // MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
     mynode->x = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.x,MAX_MORTAR_DERIVATIVES,MAX_MORTAR_DERIVATIVES/2+3*i+0),
                           MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+0);
     mynode->y = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.y,MAX_MORTAR_DERIVATIVES,MAX_MORTAR_DERIVATIVES/2+3*i+1),
                           MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+1);
     mynode->z = MadDouble(Eigen::AutoDiffScalar<DerivativeType>(node.z,MAX_MORTAR_DERIVATIVES,MAX_MORTAR_DERIVATIVES/2+3*i+2),
                           MAX_MORTAR_DERIVATIVES, MAX_MORTAR_DERIVATIVES/2+3*i+2);
#endif
     MadSlaveCoords[SlaveFace->GetNode(i)] = mynode;
   }

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
/*   // old version
     MadNormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, SlaveFacet(Facets, i), MadSlaveCoords,
                                                                             MadSlaveCoords, ngp, offset);
     MadNormalGeoGaps -= SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, MasterFacet(Facets, i), MadSlaveCoords,
                                                                             MadMasterCoords, ngp); */
     // new version
     MadNormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, SlaveFacet(Facets, i), MasterFacet(Facets, i),
                                                                             MadSlaveCoords, MadSlaveCoords, MadMasterCoords, ngp, offset);
   }

   NormalGeoGaps.reset(nMortarShapeFct,0.0);
   M.setNewSize(nMortarShapeFct,3*SlaveFace->nNodes());
   N.setNewSize(nMortarShapeFct,3*MasterFace->nNodes());
#if defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD) || defined(MORTAR_AUTO_DIFF_EIGEN_FAD_FAD)
   if(!dM) dM = new FullM[nMortarShapeFct];
   if(!dN) dN = new FullM[nMortarShapeFct];
   for(int i=0; i<nMortarShapeFct; ++i) { 
     dM[i].setNewSize(3*SlaveFace->nNodes(),MAX_MORTAR_DERIVATIVES);
     dN[i].setNewSize(3*MasterFace->nNodes(),MAX_MORTAR_DERIVATIVES);
   }
#endif
   for(int i = 0; i < nMortarShapeFct; ++i) {
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
     NormalGeoGaps[i] = MadNormalGeoGaps[i].val();
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
     Sacado::RadVec::ADvar< Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> >::Outvar_Gradcomp(MadNormalGeoGaps[i]);
     NormalGeoGaps[i] = MadNormalGeoGaps[i].val().val();
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
     NormalGeoGaps[i] = MadNormalGeoGaps[i].value();
#else // MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
     NormalGeoGaps[i] = MadNormalGeoGaps[i].value().value();
#endif
     for(int j = 0; j < SlaveFace->nNodes(); ++j) {
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
       for(int k = 0; k < 3; ++k)
         M[i][3*j+k] = MadNormalGeoGaps[i].dx(MAX_MORTAR_DERIVATIVES/2+3*j+k);
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
       M[i][3*j+0] = MadSlaveCoords[SlaveFace->GetNode(j)]->x.adj().val();
       M[i][3*j+1] = MadSlaveCoords[SlaveFace->GetNode(j)]->y.adj().val();
       M[i][3*j+2] = MadSlaveCoords[SlaveFace->GetNode(j)]->z.adj().val();
       for(int l = 0; l < MAX_MORTAR_DERIVATIVES; ++l) {
         dM[i][3*j+0][l] = MadSlaveCoords[SlaveFace->GetNode(j)]->x.adj().dx(l);
         dM[i][3*j+1][l] = MadSlaveCoords[SlaveFace->GetNode(j)]->y.adj().dx(l);
         dM[i][3*j+2][l] = MadSlaveCoords[SlaveFace->GetNode(j)]->z.adj().dx(l);
       }
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
       for(int k = 0; k < 3; ++k)
         M[i][3*j+k] = MadNormalGeoGaps[i].derivatives()[MAX_MORTAR_DERIVATIVES/2+3*j+k];
#else // MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
       for(int k = 0; k < 3; ++k) {
         M[i][3*j+k] = MadNormalGeoGaps[i].derivatives()[MAX_MORTAR_DERIVATIVES/2+3*j+k].value();
         for(int l = 0; l < MAX_MORTAR_DERIVATIVES; ++l)
           dM[i][3*j+k][l] = MadNormalGeoGaps[i].derivatives()[MAX_MORTAR_DERIVATIVES/2+3*j+k].derivatives()[l];
       }
#endif
     }

     for(int j = 0; j < MasterFace->nNodes(); ++j) {
#ifdef MORTAR_AUTO_DIFF_SACADO_FAD
       for(int k = 0; k < 3; ++k)
         N[i][3*j+k] = -MadNormalGeoGaps[i].dx(3*j+k);
#elif defined(MORTAR_AUTO_DIFF_SACADO_RAD_FAD)
       N[i][3*j+0] = -MadMasterCoords[MasterFace->GetNode(j)]->x.adj().val();
       N[i][3*j+1] = -MadMasterCoords[MasterFace->GetNode(j)]->y.adj().val();
       N[i][3*j+2] = -MadMasterCoords[MasterFace->GetNode(j)]->z.adj().val();
       for(int l = 0; l < MAX_MORTAR_DERIVATIVES; ++l) {
         dN[i][3*j+0][l] = MadMasterCoords[MasterFace->GetNode(j)]->x.adj().dx(l);
         dN[i][3*j+1][l] = MadMasterCoords[MasterFace->GetNode(j)]->y.adj().dx(l);
         dN[i][3*j+2][l] = MadMasterCoords[MasterFace->GetNode(j)]->z.adj().dx(l);
       }
#elif defined(MORTAR_AUTO_DIFF_EIGEN_FAD)
       for(int k = 0; k < 3; ++k)
         N[i][3*j+k] = -MadNormalGeoGaps[i].derivatives()[3*j+k];
#else // MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
       for(int k = 0; k < 3; ++k) {
         N[i][3*j+k] = -MadNormalGeoGaps[i].derivatives()[3*j+k].value();
         for(int l = 0; l < MAX_MORTAR_DERIVATIVES; ++l)
           dN[i][3*j+k][l] = -MadNormalGeoGaps[i].derivatives()[3*j+k].derivatives()[l];
       }
#endif
     }
   }

   for (int i=0; i<MasterFace->nNodes(); ++i)
     delete MadMasterCoords[MasterFace->GetNode(i)];
   for (int i=0; i<SlaveFace->nNodes(); ++i) 
     delete MadSlaveCoords[SlaveFace->GetNode(i)];

#ifdef MORTAR_AUTO_DIFF_SACADO_RAD_FAD
   Sacado::RadVec::ADvar< Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> >::aval_reset();
#endif
}
#endif

/*
// EXPERIMENTAL
// Alternative implementation that do not store the (triangular) facets but
// construct them on-the-fly when performing the integrations
// This may save memory as the nodes of the facets are no more duplicated and stored
// Also the FFIPolygon object may only references the ACME_FFI_Data array instead of
// extracting and storing its revelant data; in which case the ACME_FFI_Data may
// must NOT be destroyed before the FFIPolygon objects are used. This approach should
// be OK in case of multi-threads because as the ACME_FFI_Data array is only read, no 
// false-sharing should happen. 

FullM
FFIPolygon::IntegrateOnMaster_MasterShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the MASTER SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

   FullM MatShapeFctProd(nMortarShapeFct,nMasterShapeFct);
   MatShapeFctProd.zero();

   // Loop over triangularization
  //double m1M[2], m2M[2], m3M[2];
  //double m1S[2], m2S[2], m3S[2];
  int offsetm1, offsetm2;
  TriFacet MasterFacet;
  for(int i=0, ii=GetnFacets(); i<ii; i++) {
     offsetm1 = 4*i; offsetm2 = 4*((4*(i+1))%ii);

     // !! Swap Master & Slave (ACME Master <=> our Slave, ...) !!
     // !! Swap ordering to have the centroid as the first vertex of the TriFacet !! 
     MasterFacet.SetTriFacet(MasterFace, SlaveCentroid, &ACME_FFI_LocalCoordData[offsetm1], 
                                                        &ACME_FFI_LocalCoordData[offsetm2]);

     MatShapeFctProd += MasterFacet.IntegrateShapeFctProduct(MortarEl, MasterFacet, cs, ngp);
  }

  return MatShapeFctProd;
}
*/
