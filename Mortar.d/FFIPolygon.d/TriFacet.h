// --------------------------------------------------------------	
// HB - 05/29/03
// --------------------------------------------------------------	
#ifndef _TRIFACET_H_ 
#define _TRIFACET_H_

class FaceElement;
class MortarElement;
class CoordSet;

template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
#ifdef HB_NORMAL_GEOM_GAP
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
#endif 

class TriFacet {
  private:
	double Area;                     // (approximate) area 
        double LocalCoordOnFaceEl[3][2]; // coord. of the vertices in the ref. face el.
        FaceElement* FaceEl;             // ptr to the associated face el. 

  public:
        // Constructors
        // ~~~~~~~~~~~~
        TriFacet();
        TriFacet(FaceElement*, const double* m1, const double* m2, const double* m3);

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

        double GetArea() { return Area; };

        FaceElement* GetPtrFaceEl() { return FaceEl; }

        // Set methods
        // ~~~~~~~~~~~ 
        void SetTriFacet(FaceElement*, const double* m1, const double* m2, const double* m3);
        void SetArea(double);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print();
      
        // Mapping methods
        // ~~~~~~~~~~~~~~~ 
        double MappingJacobian();
        double MappingJacobian(const double* m) { return MappingJacobian(); }
        void LocalToLocalCoordOnFaceEl(const double* m, double*  mOnFaceEl);
        double GetJacobianOnFaceEl(const double* m, CoordSet& cs);
        double GetIsoParamMappingNormalAndJacobianOnFaceEl(double* Normal, const double* m, CoordSet& cs);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
        FullM IntegrateShapeFctProduct(MortarElement*, TriFacet&, CoordSet&, int ngp=6);

        FullM IntegrateNormalShapeFctProduct(MortarElement*, TriFacet&, CoordSet& , int ngp=6);

#ifdef HB_NORMAL_GEOM_GAP
        // Compute normal "geometrical" gap
        // EXPERIMENTAL
        Vector IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, TriFacet& FriendTriFacet, CoordSet &cs, int ngp=6);
#endif
};
#endif

