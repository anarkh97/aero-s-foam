// ---------------------------------------------------------------- 
// HB - 05/29/03
// HB - Modification 03/02/04
// ---------------------------------------------------------------- 
#ifndef _FFIPOLYGON_H_ 
#define _FFIPOLYGON_H_

// STL
#include <vector>
#include <utility> // std::pair

// FEM
#include <Math.d/matrix.h>
#include <Mortar.d/FFIPolygon.d/TriFacet.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

class CoordSet;
class FaceElement;
class MortarElement;

template<class Scalar>
class FFIPolygon {
  private:
	double Area;                 // (approximate) area 
        int nVertices;               // number of polygon vertices

        FaceElement* MasterFace;     // ptr to the associated master face el.
        FaceElement* SlaveFace ;     // ptr to the associated salve face el.

        typedef std::pair<TriFacet<Scalar>, TriFacet<Scalar> > Facet_pair_t; // first  -> master
                                                                             // second -> slave
        std::vector<Facet_pair_t> Facets;

        FullM M;                     // store the (FFI contibution to) M matrix (Mortar-Slave)
        FullM N;                     // store the (FFI contibution to) N matrix (Mortar-Master)

        Vector NormalGeoGaps;        // store the (FFI contibution to) the "geometrical" normal gaps (Slave-Master)

        FullM *dM;                   // store the (FFI contibution to) derivative of the M matrix (Mortar-Slave)
        FullM *dN;                   // store the (FFI contibution to) derivative of the N matrix (Mortar-Master)
#ifdef HB_ACME_FFI_DEBUG
        Scalar (*VertexLlCoordOnSFaceEl)[2];
        Scalar (*VertexLlCoordOnMFaceEl)[2];
#endif
        // Helper methods
        // ~~~~~~~~~~~~~~
        static TriFacet<Scalar>& MasterFacet(std::vector<Facet_pair_t>& FacetSet, size_t i) 
        { 
          return FacetSet[i].first; 
        }
        static TriFacet<Scalar>& SlaveFacet(std::vector<Facet_pair_t>& FacetSet, size_t i) 
        { 
          return FacetSet[i].second; 
        }

  public:
        // Constructors
        // ~~~~~~~~~~~~
        FFIPolygon();
        FFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data);
        
	// Destructor 
        // ~~~~~~~~~~
	~FFIPolygon();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Set methods
        // ~~~~~~~~~~~ 
        void SetFFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data);
        void SetPtrMasterFace(FaceElement*); 
	void SetPtrSlaveFace( FaceElement*);
	
	// Get methods  
        // ~~~~~~~~~~~
        int GetnVertices() { return nVertices; }

        int GetnFacets() { return Facets.size(); }

        FaceElement* GetPtrMasterFace() { return MasterFace; }

        FaceElement* GetPtrSlaveFace() { return SlaveFace; }

        FullM* GetPtrM() { return &M; }

        FullM* GetPtrN() { return &N; }

        Vector* GetPtrNormalGeoGaps() { return &NormalGeoGaps; }

        FullM* GetdM() { return dM; }

        FullM* GetdN() { return dN; }

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print(); 
	void PrintM();
	void PrintN();
#ifdef HB_ACME_FFI_DEBUG
	void PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId);
	void PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId);
	void PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode);
#endif
	// Triangularization methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~
	void CreateTriangularization(Scalar*);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
	FullM IntegrateOnMaster_MasterShapeFctProduct(MortarElement*, CoordSet&, int ngp=6);
	FullM IntegrateOnMaster_SlaveShapeFctProduct( MortarElement*, CoordSet&, int ngp=6);
	
        FullM IntegrateOnSlave_MasterShapeFctProduct(MortarElement*, CoordSet&, int ngp=6);
	FullM IntegrateOnSlave_SlaveShapeFctProduct( MortarElement*, CoordSet&, int ngp=6);

	FullM IntegrateOnSlave_MasterNormalShapeFctProduct(MortarElement* , CoordSet &, int ngp=6);
        FullM IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp=6);

        FullM IntegrateOnSlave_MasterGradNormalShapeFctProduct(MortarElement*, CoordSet&, CoordSet&, int ngp=6);
        FullM IntegrateOnSlave_SlaveGradNormalShapeFctProduct(MortarElement*, CoordSet&, CoordSet&, int ngp=6);

        void ComputeM(MortarElement*, CoordSet&, int ngp=6);
        void ComputeN(MortarElement*, CoordSet&, int ngp=6);

        void ComputeNormalM(MortarElement*, CoordSet&, int ngp=6);
        void ComputeNormalN(MortarElement*, CoordSet&, int ngp=6);

        void ComputeGradNormalM(MortarElement*, CoordSet&, CoordSet&, int ngp=6);
        void ComputeGradNormalN(MortarElement*, CoordSet&, CoordSet&, int ngp=6);

	void ComputeNormalGeoGap(MortarElement* MortarEl, CoordSet&, CoordSet&, int ngp=6, double offset=0.);

        // Space/memory allocation/desallocation methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
};

#ifdef _TEMPLATE_FIX_
  #include <Mortar.d/FFIPolygon.d/FFIPolygon.C>
#endif

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
class FFIPolygon {
  private:
        double Area;                 // (approximate) area
        int nVertices;               // number of polygon vertices

        FaceElement* MasterFace;     // ptr to the associated master face el.
        FaceElement* SlaveFace ;     // ptr to the associated salve face el.

        double MasterCentroid[2];
        double SlaveCentroid[2];

        double* ACME_FFI_Data;

        FullM M;                     // store the (FFI contibution to) M matrix (Mortar-Slave)
        FullM N;                     // store the (FFI contibution to) N matrix (Mortar-Master)

        Vector NormalGeoGaps;        // store the (FFI contibution to) the "geometrical" normal gaps (Slave-Master)

#ifdef HB_ACME_FFI_DEBUG
        double (*VertexLlCoordOnSFaceEl)[2];
        double (*VertexLlCoordOnMFaceEl)[2];
#endif
        // Helper methods
        // ~~~~~~~~~~~~~~
        double* ACME_FFI_LocalCoordData() { return &ACME_FFI_Data[2*nVertices+1]; }

  public:
        // Constructors
        // ~~~~~~~~~~~~
        FFIPolygon();
        FFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data);

        // Destructor
        // ~~~~~~~~~~
        ~FFIPolygon();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Set methods
        // ~~~~~~~~~~~
        void SetFFIPolygon(FaceElement* MasterFaceEl, FaceElement* SlaveFaceEl, int nVert, double* ACME_FFI_Data);
        void SetPtrMasterFace(FaceElement*);
        void SetPtrSlaveFace( FaceElement*);

        // Get methods
        // ~~~~~~~~~~~
        int GetnVertices() { return nVertices; }

        int GetnFacets() { return nVertices; }

        FaceElement* GetPtrMasterFace() { return MasterFace; }

        FaceElement* GetPtrSlaveFace() { return SlaveFace; }

        FullM* GetPtrM() { return &M; }

        FullM* GetPtrN() { return &N; }

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print();
        void PrintM();
        void PrintN();
#ifdef HB_ACME_FFI_DEBUG
        void PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId);
        void PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId);
        void PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode);
#endif
        // Triangularization methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        void CreateTriangularization(double*);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~
        FullM IntegrateOnMaster_MasterShapeFctProduct(MortarElement*, CoordSet&, int ngp=6);
        FullM IntegrateOnMaster_SlaveShapeFctProduct( MortarElement*, CoordSet&, int ngp=6);

        FullM IntegrateOnSlave_MasterShapeFctProduct(MortarElement*, CoordSet&, int ngp=6);
        FullM IntegrateOnSlave_SlaveShapeFctProduct( MortarElement*, CoordSet&, int ngp=6);

        FullM IntegrateOnSlave_MasterNormalShapeFctProduct(MortarElement* , CoordSet &, int ngp=6);
        FullM IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarElement* MortarEl, CoordSet &cs, int ngp=6);

        void ComputeM(MortarElement*, CoordSet&, int ngp=6);
        void ComputeN(MortarElement*, CoordSet&, int ngp=6);

        void ComputeNormalM(MortarElement*, CoordSet&, int ngp=6);
        void ComputeNormalN(MortarElement*, CoordSet&, int ngp=6);

        void ComputeNormalGeoGap(MortarElement* MortarEl, CoordSet &cs, int ngp=6);

        // Space/memory allocation/desallocation methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
};
#endif

*/
