// ---------------------------------------------------------------- 
// HB - 05/06/03
// ---------------------------------------------------------------- 
#ifndef _FACEELEMENT_H_
#define _FACEELEMENT_H_

#include <iostream>

// STL 
#include <map>

// FEM headers
#include <Element.d/Element.h>
#include <Utils.d/resize_array.h>
#include <Math.d/matrix.h>

//class CoordSet;
class FFIPolygon;
class DofSetArray;
//class State;
//class GeomState;
struct InterpPoint;
class Connectivity;

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

class FaceElement {

	double Area; 

	//int nFFI;	
	//ResizeArray<FFIPolygon*> FFI; 

        //int ElemId;
        //int PtrElem;
  public:
        // public data 
	// ~~~~~~~~~~~	
	//enum FaceElemType{QUADFACEL4=1, QUADFACEQ8, TRIFACEL3, TRIFACEQ6};
        enum {QUADFACEL4=1, QUADFACEQ8, TRIFACEL3, TRIFACEQ6, SHELLQUADFACEL4, SHELLTRIFACEL3, QUADFACEQ9, QUADFACEC12, TRIFACEC10, POINTFACE};
      
	// Constructors 
	// ~~~~~~~~~~~~
	FaceElement():Area(0.0) { }
	//FaceElement():FFI(0) { Area = 0.0; nFFI = 0; }
	//FaceElement() { Area = 0.0; nFFI = 0; }
	//FaceElement();
        
        // for future
        //virtual FaceElement* clone(); 

        // Setup & Update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        virtual void Renumber(std::map<int,int>& OldToNewNodeIds)=0;
 	
        // Get methods
        // ~~~~~~~~~~~
	// -> pure interface methods
	virtual int nNodes()=0;    
        virtual void GetNodes(int*, int* renumTable=0)=0;
        virtual void GetNodes(int*, std::map<int,int>& renumTable)=0;
        virtual int GetNode(int)=0;

        virtual int GetNodeIndex(int)=0;
	virtual int GetFaceElemType()=0;
        
        // -> for the case of quadratic elements which are NOT 
        //    yet supported by ACME for FFI
        //    basicaly, these methods will return the "linear" nodes
        //    of the element: for example, in the case of a FaceTri6 element
        //    these methods will return the 3 first nodes which are in fact
        //    the 3 vertices nodes of the underlaying FaceTri3 element
        // -> So we will perform the FFI search using the underlaying "linear"
        //    face element: FaceQuad8 -> FaceQuad4 
        //                  FaceTri6  -> FaceTri3
        //    THIS IS AN APPROXIMATION !! IT IS OK FOR STRAIGHT EDGES FACE 
        //    ELEMENT, BUT WE MAKE ERRORS IN THE CASE OF CURVED EDGES.
        //    NOTE THAT THIS IS ONLY VALID FOR THE "GEOMETRIC" FFI SEARCH.
        //    FOR THE MORTAR INTEGRATION, WE WILL USE THE "TRUE" (QUADRATIC)
        //    MAPPING.
        // -> in the case of "linear" face element, these methods are equivalent
        //    to nNodes(), GetNode(int) & GetNodes(int*).
        virtual int nVertices()=0;
        virtual int GetVertex(int)=0;
        virtual void GetVertices(int*, int* renumTable=0)=0;
        virtual void GetVertices(int*, std::map<int,int>& renumTable)=0;

	// -> for ACME
#ifdef USE_ACME
	virtual ContactSearch::ContactFace_Type GetACMEFaceElemType()=0;
#else
	virtual int GetACMEFaceElemType()=0;
#endif        
        // -> see above, the discussion about the "quadratic" face elements.	
        //    This method will return the ACME face element type used for 
        //    the FFI search: FaceQuad8 -> FaceQuad4
        //                    FaceTri6  -> FaceTri3
        //    For "linear" face elements, this method is equivalent to
        //    GetACMEFaceElemType(). 
#ifdef USE_ACME
        virtual ContactSearch::ContactFace_Type GetACMEFFIFaceElemType()=0;
#else
        virtual int GetACMEFFIFaceElemType()=0;
#endif

        // Mapping & shape fct methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// -> pure interface methods
 	//virtual void GetShapeFct(double* Shape, double* m)=0;
 	virtual void GetShapeFctVal(double* Shape, double* m)=0;

        virtual double GetJacobian(double* m, CoordSet &cs)=0;
        virtual double GetIsoParamMappingNormalAndJacobian(double*, double*, CoordSet &cs)=0;

        virtual void LocalToGlobalCoord(double*, double*, CoordSet&)=0;

        virtual FullM ScalarMass(CoordSet& , double rho, int ngp)=0;
        virtual void IntegrateShapeFcts(double*, CoordSet&, double rho, int ngp)=0;

        virtual double* ViewRefCoords();
        virtual void GetdNormal(double dNormal[][3], double* m, CoordSet& cs) { 
          std::cerr << "FaceElement::GetdNormal not implemented\n"; exit(-1); }

	// Print, display ... methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~
	virtual void print()=0;

        virtual int numDofs() {fprintf(stderr,"function numDofs() undefined for this type of element!\n"); return 0;}	
        virtual int* dofs(DofSetArray &, int*, int*) {fprintf(stderr,"function dofs(...) undefined for this type of element!\n"); return 0;}
        virtual void computeDisp(CoordSet&, State&, const InterpPoint&, double*, GeomState*, int*) {
          fprintf(stderr,"function computeDisp(...) undefined for this type of element!\n");}
        virtual void getFlLoad(const InterpPoint&, double*, double*) {
          fprintf(stderr,"function computeDisp(...) undefined for this type of element!\n");}

	// FFI methods
	// ~~~~~~~~~~~
	//void AddPtrFFI(FFIPolygon*);
	//int nFFIs();
	//FFIPolygon* GetPtrFFI(int);
	//void printFFI();

        int findEle(Connectivity *nodeToElem, int *eleTouch,
                    int *eleCount, int myNum, int *fnId);

};
#endif

