// --------------------------------------------------------------	
// HB - 05/29/03
// --------------------------------------------------------------	
#ifndef _TRIFACET_H_ 
#define _TRIFACET_H_

class FaceElement;
class MortarElement;

template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

template <class Scalar>
class TriFacet {
  private:
	//double Area;                     // (approximate) area 
        Scalar LocalCoordOnFaceEl[3][2]; // coord. of the vertices in the ref. face el.
        FaceElement* FaceEl;             // ptr to the associated face el. 

  public:
        // Constructors
        // ~~~~~~~~~~~~
        TriFacet();
        TriFacet(FaceElement*, const Scalar* m1, const Scalar* m2, const Scalar* m3);

        // Destructor
        // ~~~~~~~~~~
        ~TriFacet();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Get methods  
        // ~~~~~~~~~~~~
        int nVertices() { return 3; }

        int nNodes() { return nVertices(); }

        //double GetArea() { return Area; };

        FaceElement* GetPtrFaceEl() { return FaceEl; }

        // Set methods
        // ~~~~~~~~~~~ 
        void SetTriFacet(FaceElement*, const Scalar* m1, const Scalar* m2, const Scalar* m3);
        //void SetArea(double);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print();
      
        // Mapping methods
        // ~~~~~~~~~~~~~~~ 
        Scalar MappingJacobian();
        Scalar MappingJacobian(const Scalar* m) { return MappingJacobian(); }
        void LocalToLocalCoordOnFaceEl(const double* m, Scalar* mOnFaceEl);
        template<typename CoordSetT>
        Scalar GetJacobianOnFaceEl(const Scalar* m, CoordSetT& cs);
        template<typename CoordSetT>
        Scalar GetIsoParamMappingNormalAndJacobianOnFaceEl(Scalar* Normal, const Scalar* m, CoordSetT& cs);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
        FullM IntegrateShapeFctProduct(MortarElement*, TriFacet&, CoordSet&, int ngp=6);
        FullM IntegrateNormalShapeFctProduct(MortarElement*, TriFacet&, CoordSet&, int ngp=6);
        FullM IntegrateGradNormalShapeFctProduct(MortarElement*, TriFacet&, CoordSet&, CoordSet&,
                                                 TriFacet&, CoordSet&, int ngp=6);

        // Compute normal "geometrical" gap
/* old version
        template<typename CoordSetT>
        GenVector<Scalar> IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, TriFacet& FriendTriFacet,
                                                        CoordSetT& cs, CoordSetT& cs1, int ngp=6, double offset=0.); */
        template<typename CoordSetT>
        GenVector<Scalar> IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, TriFacet& FriendFacetB, TriFacet& FriendFacetC,
                                                        CoordSetT& cs, CoordSetT& csB, CoordSetT& csC, int ngp=6,
                                                        double offset=0.);
};

#ifdef _TEMPLATE_FIX_
  #include <Mortar.d/FFIPolygon.d/TriFacet.C>
#endif

#endif
