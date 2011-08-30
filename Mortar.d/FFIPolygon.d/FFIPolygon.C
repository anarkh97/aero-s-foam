// Std C/C++ headers
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/FFIPolygon.d/TriFacet.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
FFIPolygon::FFIPolygon()
: Facets()
{
  Initialize();
}

FFIPolygon::FFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data)
: Facets()
{
  Initialize();
  SetFFIPolygon(MasterFaceEl, SlaveFaceEl, nVert, ACME_FFI_Data);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
FFIPolygon::~FFIPolygon()
{
#ifdef HB_ACME_FFI_DEBUG
  if(VertexLlCoordOnSFaceEl) { delete VertexLlCoordOnSFaceEl; }
  if(VertexLlCoordOnMFaceEl) { delete VertexLlCoordOnMFaceEl; }
#endif
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
void 
FFIPolygon::Initialize()
{
  Area      = 0.0;
  nVertices = 0;
  MasterFace= 0;
  SlaveFace = 0;

#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = 0;
  VertexLlCoordOnMFaceEl = 0;
#endif
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
void
FFIPolygon::SetPtrMasterFace(FaceElement* PtrMasterFace) { MasterFace = PtrMasterFace; }

void
FFIPolygon::SetPtrSlaveFace(FaceElement* PtrSlaveFace) { SlaveFace = PtrSlaveFace; }

void
FFIPolygon::SetFFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, 
                          int nVert, double* ACME_FFI_Data)
{
  nVertices  = nVert;
  MasterFace = MasterFaceEl;
  SlaveFace  = SlaveFaceEl;

  // Create & fill polygon triangularization
  int LocalCoordOffset = 2*nVertices+1;
  //fprintf(stderr," In FFIPolygon::SetFFIPolygon, LocalCoordOffset = %d\n",LocalCoordOffset);
 
  double* ACME_FFI_LocalCoordData = &ACME_FFI_Data[LocalCoordOffset];
  CreateTriangularization(ACME_FFI_LocalCoordData);
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
void
FFIPolygon::Print()
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

void
FFIPolygon::PrintM()
{
  GetPtrSlaveFace()->print(); 
  GetPtrSlaveFace()->print(); 
  M.print("M[FFI] = ");
}

void
FFIPolygon::PrintN()
{
  GetPtrSlaveFace()->print(); 
  GetPtrMasterFace()->print(); 
  N.print("N[FFI] = ");
}

#ifdef HB_ACME_FFI_DEBUG
void 
FFIPolygon::PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnSFaceEl[iVert][0],VertexLlCoordOnSFaceEl[iVert][1]};
    SlaveFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

void
FFIPolygon::PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnMFaceEl[iVert][0],VertexLlCoordOnMFaceEl[iVert][1]};
    MasterFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

void
FFIPolygon::PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode)
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
void
FFIPolygon::CreateTriangularization(double* ACME_FFI_LocalCoordData)
{
  //fprintf(stderr," In FFIPolygon::CreateTriangularization:\n");
  
  // Allocate space for Facets
  // By default assume as many triangular facet as polygon vertices
  // -> can be overwritten in for optimization
  Facets.assign(nVertices, Facet_pair_t(TriFacet(), TriFacet()));

  // Create master & slave polygon centroids
  double MasterCentroid[2] = {0.0, 0.0};
  double SlaveCentroid[2]  = {0.0, 0.0};
#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = new double[nVertices][2];
  VertexLlCoordOnMFaceEl = new double[nVertices][2];
#endif
  for(int iVert=0; iVert<nVertices; ++iVert) {
     int offset = 4*iVert;
     //fprintf(stderr," --------------------------------------- \n");
     //fprintf(stderr," ACME_FFI_LocalCoordData[%2d] = %e\n",offset  ,ACME_FFI_LocalCoordData[offset  ]);
     //fprintf(stderr," ACME_FFI_LocalCoordData[%2d] = %e\n",offset+1,ACME_FFI_LocalCoordData[offset+1]);
     //fprintf(stderr," ACME_FFI_LocalCoordData[%2d] = %e\n",offset+2,ACME_FFI_LocalCoordData[offset+2]);
     //fprintf(stderr," ACME_FFI_LocalCoordData[%2d] = %e\n",offset+3,ACME_FFI_LocalCoordData[offset+3]);
#ifdef HB_ACME_FFI_DEBUG
     // PJSA swap master and slave (ACME Master <=> our Slave, ...)
     VertexLlCoordOnMFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset  ];
     VertexLlCoordOnMFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+1];
     VertexLlCoordOnSFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset+2];
     VertexLlCoordOnSFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+3];
#endif 
     //SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset  ]; 
     //SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+1]; 
     //MasterCentroid[0] += ACME_FFI_LocalCoordData[offset+2]; 
     //MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+3]; 
     //HB: Swap Master & Slave (ACME Master <=> our Slave, ...) 
     //    do it here so that the code below is more "natural"
     MasterCentroid[0] += ACME_FFI_LocalCoordData[offset  ]; 
     MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+1]; 
     SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset+2]; 
     SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+3]; 
  }
  SlaveCentroid[0]  /= nVertices; SlaveCentroid[1]  /= nVertices;
  MasterCentroid[0] /= nVertices; MasterCentroid[1] /= nVertices;

  //fprintf(stderr," SlaveCentroid : x = %e, y = %e\n",SlaveCentroid[0], SlaveCentroid[1]);
  //fprintf(stderr," MasterCentroid: x = %e, y = %e\n",MasterCentroid[0], MasterCentroid[1]);

  // Create master & slave triangularization
  for(size_t i=0, ii=GetnFacets(), offsetm1=0; i<ii; i++, offsetm1+=4) {
     size_t offsetm2 = (i==ii-1) ? 0 : offsetm1+4;
     //HB: note that ACME Master <=> our Slave, ... 
     MasterFacet(Facets, i).SetTriFacet(MasterFace, MasterCentroid , &ACME_FFI_LocalCoordData[offsetm1],
                                                                     &ACME_FFI_LocalCoordData[offsetm2]);
     SlaveFacet(Facets , i).SetTriFacet(SlaveFace , SlaveCentroid,   &ACME_FFI_LocalCoordData[offsetm1+2],
                                                                     &ACME_FFI_LocalCoordData[offsetm2+2]);
  } 
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FFIPolygon::ComputeArea()
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
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     MatShapeFctProd += MasterFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                        cs, ngp);
   }

   return MatShapeFctProd;
}

FullM
FFIPolygon::IntegrateOnMaster_SlaveShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
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

FullM
FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
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

   //cerr << "In FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct" << endl;
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

FullM
FFIPolygon::IntegrateOnSlave_SlaveShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
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
   
   //cerr << "In FFIPolygon::IntegrateOnSlave_SlaveShapeFctProduct" << endl;
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

void 
FFIPolygon::ComputeM(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Slave shape fcts product mass-like matrix M
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon::IntegrateOnSlave_SlaveShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
//         -> see FFIPolygon::IntegrateOnSlave_SlaveShapeFctProduct for more details
// *****************************************************************************************************
{
  M = IntegrateOnSlave_SlaveShapeFctProduct(MortarEl, cs, ngp);
  //M.print();
}

void 
FFIPolygon::ComputeN(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterShapeFctProduct(MortarEl, cs, ngp);
  //N.print();
}

FullM
FFIPolygon::IntegrateOnSlave_MasterNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
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

FullM
FFIPolygon::IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp)
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

void
FFIPolygon::ComputeNormalM(MortarElement* MortarEl, CoordSet &cs, int ngp)
{
  M = IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarEl, cs, ngp);
}

void
FFIPolygon::ComputeNormalN(MortarElement* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygon::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterNormalShapeFctProduct(MortarEl, cs, ngp);
}

FullM
FFIPolygon::IntegrateOnSlave_MasterGradNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
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

FullM
FFIPolygon::IntegrateOnSlave_SlaveGradNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
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

void
FFIPolygon::ComputeGradNormalM(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
  M = IntegrateOnSlave_SlaveGradNormalShapeFctProduct(MortarEl, SlaveCoords, MasterCoords, ngp);
}

void
FFIPolygon::ComputeGradNormalN(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
  N = IntegrateOnSlave_MasterGradNormalShapeFctProduct(MortarEl, SlaveCoords, MasterCoords, ngp);
}

void
FFIPolygon::ComputeNormalGeoGap(MortarElement* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp)
{
   int nMortarShapeFct = MortarEl->nNodes();

   NormalGeoGaps.reset(nMortarShapeFct,0.0);

   // Loop over triangularization
   for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
     NormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, SlaveFacet(Facets, i), SlaveCoords,
                                                                          SlaveCoords, ngp);
     NormalGeoGaps -= SlaveFacet(Facets, i).IntegrateNormalGeoGagsProduct(MortarEl, MasterFacet(Facets, i), SlaveCoords,
                                                                          MasterCoords, ngp);
   }
   //NormalGeoGaps.print(" NormalGeoGaps[FFIPolygon]");
}

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
