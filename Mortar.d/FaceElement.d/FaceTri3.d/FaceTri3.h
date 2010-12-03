// ----------------------------------------------------------------
// HB - 08/15/03
// ----------------------------------------------------------------
#ifndef _FACETRI3_H_
#define _FACETRI3_H_

// STL
#include <map>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
//#include <Utils.d/dofset.h>
//#include <Hetero.d/FlExchange.h>
//#include <Element.d/State.h>
class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

class FaceTri3: public FaceElement {
  private:
        int Nodes[3];
        static double RefCoords[3][2]; // coords of the nodes in the ref./parametric domain

  public:
        // Constructors
	// ~~~~~~~~~~~~
        FaceTri3(int *);

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual methods
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
	void   GetShapeFct(double* Shape, double* m);
        double GetShapeFctAndJacobian(double* Shape, double* m, CoordSet& cs);
	double GetJacobian(CoordSet& cs);
	
	// -> implementation of virtual fcts
        double* ViewRefCoords();

	// -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double* M, double* m, CoordSet& cs);
        void   GetShapeFctVal(double* Shape, double* m);
	double GetJacobian(double* m, CoordSet& cs);
        double GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs);

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

        int numDofs() {return 9;}
        int* dofs(DofSetArray &dsa, int *p, int *fnId); 
        void computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, GeomState*, int *fnId); 
        void getFlLoad(CoordSet&, const InterpPoint &ip, double *flF, double *resF, GeomState*); 

};

#endif
