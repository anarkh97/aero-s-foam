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
        template<typename Scalar>
	  void   GetShapeFct(Scalar* Shape, Scalar* m);
        template<typename Scalar>
          void   GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar, typename CoordSetT>
          Scalar GetJacobian(CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          Scalar GetShapeFctAndJacobian(Scalar* Shape, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   LocalToGlobalCoordTemp(Scalar* M, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetIsoParamMappingNormalJacobianProductTemp(Scalar* JNormal, Scalar* m, CoordSetT& cs);
	
	// -> implementation of virtual fcts
        double* ViewRefCoords();
        void GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs);

	// -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double* M, double* m, CoordSet& cs);
        void   GetShapeFctVal(double* Shape, double* m);
	double GetJacobian(double* m, CoordSet& cs);
        double GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs);
        void   GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs);
#if (MAX_MORTAR_DERIVATIVES > 0)
        void   GetShapeFctVal(ActiveDouble*, ActiveDouble*);
        void   LocalToGlobalCoord(ActiveDouble* M, ActiveDouble* m, MadCoordSet& cs);
        void   GetIsoParamMappingNormalJacobianProduct(ActiveDouble* JNormal, ActiveDouble* m, MadCoordSet& cs);
#endif

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
        void getFlLoad(const InterpPoint &ip, double *flF, double *resF); 

};


template<typename Scalar>
void
FaceTri3::GetShapeFct(Scalar *Shape, Scalar *m)
{
  Scalar r = m[0];
  Scalar s = m[1];
  Scalar t = 1.-r-s;

  // !! idem ACME !! 
  Shape[0] = r;
  Shape[1] = s;
  Shape[2] = t;
}

template<typename Scalar>
void
FaceTri3::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  Scalar x = m[0];
  Scalar y = m[1];

  dShapex[0] = 1;
  dShapex[1] = 0;
  dShapex[2] = -1;

  dShapey[0] = 0;
  dShapey[1] = 1;
  dShapey[2] = -1;
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri3::GetJacobian(CoordSetT &cs)
{
  // J = 2*Area = ||12 x 13||
  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  Scalar V12[3], V13[3];
  V12[0] = X[1]-X[0]; V12[1] = Y[1]-Y[0]; V12[2] = Z[1]-Z[0];
  V13[0] = X[2]-X[0]; V13[1] = Y[2]-Y[0]; V13[2] = Z[2]-Z[0];

  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  return(sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri3::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  Scalar r = m[0];
  Scalar s = m[1];
  Scalar t = 1.-r-s;

  // !! idem ACME !! 
  Shape[0] = r;
  Shape[1] = s;
  Shape[2] = t;

  return(GetJacobian(cs));
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[3], dShapey[3];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = dShapex[0]*X[0] + dShapex[1]*X[1] + dShapex[2]*X[2];
  dMdx[1] = dShapex[0]*Y[0] + dShapex[1]*Y[1] + dShapex[2]*Y[2];
  dMdx[2] = dShapex[0]*Z[0] + dShapex[1]*Z[1] + dShapex[2]*Z[2];

  dMdy[0] = dShapey[0]*X[0] + dShapey[1]*X[1] + dShapey[2]*X[2];
  dMdy[1] = dShapey[0]*Y[0] + dShapey[1]*Y[1] + dShapey[2]*Y[2];
  dMdy[2] = dShapey[0]*Z[0] + dShapey[1]*Z[1] + dShapey[2]*Z[2];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::LocalToGlobalCoordTemp(Scalar *M, Scalar *m, CoordSetT &cs)
{
  Scalar r = m[0];
  Scalar s = m[1];
  Scalar t = 1.-r-s;

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // !! idem ACME !!
  M[0] = r*X[0]+s*X[1]+t*X[2];
  M[1] = r*Y[0]+s*Y[1]+t*Y[2];
  M[2] = r*Z[0]+s*Z[1]+t*Z[2];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetIsoParamMappingNormalJacobianProductTemp(Scalar *JNormal, Scalar *m, CoordSetT &cs)
{
  // J = 2*Area = ||12 x 13||
  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  Scalar V12[3], V13[3];
  V12[0] = X[1]-X[0]; V12[1] = Y[1]-Y[0]; V12[2] = Z[1]-Z[0];
  V13[0] = X[2]-X[0]; V13[1] = Y[2]-Y[0]; V13[2] = Z[2]-Z[0];

  JNormal[0] = V12[1]*V13[2] - V12[2]*V13[1];
  JNormal[1] = V12[2]*V13[0] - V12[0]*V13[2];
  JNormal[2] = V12[0]*V13[1] - V12[1]*V13[0];
}

#endif
