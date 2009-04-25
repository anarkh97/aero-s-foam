// ----------------------------------------------------------------
// HB - 08/26/03
// ----------------------------------------------------------------
#ifndef _FACEQUAD8_H_
#define _FACEQUAD8_H_

// STL
#include <map>

#include <Mortar.d/FaceElement.d/FaceElement.h>

#ifdef USE_ACME
// ACME headers
#include "ContactSearch.h"
#endif

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

class FaceQuad8: public FaceElement {
  private:
        int Nodes[8];
        static double RefCoords[8][2]; // coords of the nodes in the ref./parametric domain

  public:
        // Constructors
	// ~~~~~~~~~~~~
        FaceQuad8(int *);

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual methods 
        void Renumber(std::map<int,int>& OldToNewNodeIds);

        // Get methods
	// ~~~~~~~~~~~
        // -> local methods 
	int  nQuad4Nodes(); 
	int  GetQuad4Node(int); 
	void GetQuad4Nodes(int*, int* renumTable=0);
	void GetQuad4Nodes(int*, std::map<int,int>& renumTable);

        // -> implementation of pure virtual methods 
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
	void   GetdShapeFct(double *dShapex, double *dShapey, double *m); 
        double GetShapeFctAndJacobian(double *, double *, CoordSet&);
        double GetIsoParamMappingNormalAndJacobian(double*, double*, CoordSet&);
        void   ComputedMdxAnddMdy(double *dMdx, double *dMdy, double *m, CoordSet &cs);
	// -> implementation of virtual fcts
        double* ViewRefCoords();

	// -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double*, double*, CoordSet&);
        void   GetShapeFctVal(double*, double*);
	double GetJacobian(double*, CoordSet&);

	// Miscelleaneous methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of virtual fcts
	//double ComputeArea(CoordSet&, const int ngp=3)

        // Mass matrix methods
	// ~~~~~~~~~~~~~~~~~~~
	// -> implementation of pure virtual fcts
	FullM ScalarMass(CoordSet& , double rho=1.0, int ngp=3);
	void  IntegrateShapeFcts(double*, CoordSet&, double rho=1.0, int ngp=3);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> local fcts
	void printNodes();

	// -> implementation of pure virtual fcts
        void print();

};
#endif
