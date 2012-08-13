#ifndef _FACEPOINT1_H_
#define _FACEPOINT1_H_

// STL
#include <map>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// FEM headers
#include <Element.d/Element.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

class FacePoint1: public FaceElement {
  private:
        int Nodes[1];
        static double RefCoords[1][2]; // coords of the nodes in the ref./parametric domain

  public:
        // Constructors
	// ~~~~~~~~~~~~
        FacePoint1(int *);

        // for future: copy constructor
        //FacePoint1(const FacePoint1 &FQ4);
        
	// Copy & clone methods
        // ~~~~~~~~~~~~~~~~~~~~
	// FaceElement* clone();

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
	void Renumber(std::map<int,int>& OldToNewNodeIds);

        // Get methods
	// ~~~~~~~~~~~
        // -> implementation of virtual fcts
	int  nNodes();
	void GetNodes(int*, int* renumTable);
	void GetNodes(int*, std::map<int,int>& renumTable);
	int  GetNode(int);
        int  GetNodeIndex(int);

	int GetFaceElemType();
#ifdef USE_ACME
        ContactSearch::ContactFace_Type GetACMEFaceElemType();
#else
	int GetACMEFaceElemType();
#endif
        // -> pure virtual method for dealing with quadratic face element
        //    (see FaceElement.h for more details)
        int  nVertices();
        int  GetVertex(int);
        void GetVertices(int*, int* renumTable);
        void GetVertices(int*, std::map<int,int>& renumTable);
    
#ifdef USE_ACME
        ContactSearch::ContactFace_Type GetACMEFFIFaceElemType();
#else
	int GetACMEFFIFaceElemType();
#endif
	// Mapping & shape fct methods        
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // -> local methods 
	void   GetShapeFct(double*, double*);
        void   GetdShapeFct(double* dShapex, double* dShapey, double* m);
	double GetShapeFctAndJacobian(double* Shape, double* m, CoordSet&);
	void   ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs);
	
	// -> implementation of virtual fcts
        double* ViewRefCoords();

	// -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double* M, double* m, CoordSet& cs);
        void   GetShapeFctVal(double*, double*);
	double GetJacobian(double*, CoordSet&);
        double GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs);
        void   GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs);

	// Miscelleaneous methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of virtual fcts
	//double ComputeArea(CoordSet&, const int ngp=2)

        // Mass matrix methods
	// ~~~~~~~~~~~~~~~~~~~
	// -> implementation of pure virtual fcts
	FullM ScalarMass(CoordSet& , double rho=1.0, int ngp=2);
	void IntegrateShapeFcts(double*, CoordSet&, double rho=1.0, int ngp=2);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> local fcts
	void printNodes();

	// -> implementation of pure virtual fcts
        void print();
};

#endif
