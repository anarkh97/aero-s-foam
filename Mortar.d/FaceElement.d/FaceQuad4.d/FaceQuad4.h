// ----------------------------------------------------------------
// HB - 05/06/03
// ----------------------------------------------------------------
#ifndef _FACEQUAD4_H_
#define _FACEQUAD4_H_

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

class FaceQuad4: public FaceElement {
  private:
        int Nodes[4];
        static double RefCoords[4][2]; // coords of the nodes in the ref./parametric domain

  public:
        // Constructors
	// ~~~~~~~~~~~~
        FaceQuad4(int *);

        // for future: copy constructor
        //FaceQuad4(const FaceQuad4 &FQ4);
        
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
        template<typename Scalar>
	  void   GetShapeFct(Scalar*, Scalar*);
        template<typename Scalar>
          void   GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar, typename CoordSetT>
	  Scalar GetShapeFctAndJacobian(Scalar* Shape, Scalar* m, CoordSetT&);
        template<typename Scalar, typename CoordSetT>
	  void   ComputedMdxAnddMdy(Scalar* dMdx, Scalar* dMdy, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   LocalToGlobalCoordTemp(Scalar* M, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetIsoParamMappingNormalJacobianProductTemp(Scalar* JNormal, Scalar* m, CoordSetT& cs);
	
	// -> implementation of virtual fcts
        double* ViewRefCoords();
        void GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs);

	// -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double* M, double* m, CoordSet& cs);
        void   GetShapeFctVal(double*, double*);
	double GetJacobian(double*, CoordSet&);
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

};

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------

template<typename Scalar>
void
FaceQuad4::GetShapeFct(Scalar *Shape, Scalar *m)
{
  Scalar x = m[0];
  Scalar y = m[1];

  Scalar d1 = 0.5*(1.0+x);
  Scalar d2 = 0.5*(1.0+y);
  Scalar d3 = 1.0-d1;
  Scalar d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;
}

template<typename Scalar>
void
FaceQuad4::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  Scalar x = m[0];
  Scalar y = m[1];
  Scalar onequart = 1./4.;

  Scalar xm = 1.-x;
  Scalar xp = 1.+x;
  Scalar ym = 1.-y;
  Scalar yp = 1.+y;

  dShapex[0] = -onequart*ym;
  dShapex[1] =  onequart*ym;
  dShapex[2] =  onequart*yp;
  dShapex[3] = -onequart*yp;

  dShapey[0] = -onequart*xm;
  dShapey[1] = -onequart*xp;
  dShapey[2] =  onequart*xp;
  dShapey[3] =  onequart*xm;
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceQuad4::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  Scalar x = m[0];
  Scalar y = m[1];

  Scalar d1 = 0.5*(1.0+x);
  Scalar d2 = 0.5*(1.0+y);
  Scalar d3 = 1.0-d1;
  Scalar d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;

  double X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  Scalar a[4], b[4], c[4];
  a[0] = (Y[1]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[1]-Z[0]);
  a[1] = (Y[1]-Y[0])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[1]-Z[0]);
  a[2] = (Y[2]-Y[3])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[2]-Z[3]);
  a[3] = (Y[2]-Y[3])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[3]);

  b[0] = (Z[1]-Z[0])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[1]-X[0]);
  b[1] = (Z[1]-Z[0])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[1]-X[0]);
  b[2] = (Z[2]-Z[3])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[2]-X[3]);
  b[3] = (Z[2]-Z[3])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[2]-X[3]);

  c[0] = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]);
  c[1] = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]);
  c[2] = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]);
  c[3] = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]);

  Scalar N[3];
  N[0] = Shape[0]*a[0]+Shape[1]*a[1]+Shape[2]*a[2]+Shape[3]*a[3];
  N[1] = Shape[0]*b[0]+Shape[1]*b[1]+Shape[2]*b[2]+Shape[3]*b[3];
  N[2] = Shape[0]*c[0]+Shape[1]*c[1]+Shape[2]*c[2]+Shape[3]*c[3];

  return(0.25*sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute dM/dx & dM/dy
  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = dShapex[0]*X[0] + dShapex[1]*X[1] + dShapex[2]*X[2] + dShapex[3]*X[3];
  dMdx[1] = dShapex[0]*Y[0] + dShapex[1]*Y[1] + dShapex[2]*Y[2] + dShapex[3]*Y[3];
  dMdx[2] = dShapex[0]*Z[0] + dShapex[1]*Z[1] + dShapex[2]*Z[2] + dShapex[3]*Z[3];

  dMdy[0] = dShapey[0]*X[0] + dShapey[1]*X[1] + dShapey[2]*X[2] + dShapey[3]*X[3];
  dMdy[1] = dShapey[0]*Y[0] + dShapey[1]*Y[1] + dShapey[2]*Y[2] + dShapey[3]*Y[3];
  dMdy[2] = dShapey[0]*Z[0] + dShapey[1]*Z[1] + dShapey[2]*Z[2] + dShapey[3]*Z[3];
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::LocalToGlobalCoordTemp(Scalar* M, Scalar* m, CoordSetT &cs)
{
  // input  : m, the local coordinates of some point P on the face
  //         cs, the global coordinates of the vertices of the face 
  // output : M, the global x,y,z coordinates of P
  Scalar Shape[4];
  GetShapeFct(Shape,m);

  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y; 
    Z[i] = cs[Nodes[i]]->z;
  }

  M[0] = Shape[0]*X[0]+Shape[1]*X[1]+Shape[2]*X[2]+Shape[3]*X[3];
  M[1] = Shape[0]*Y[0]+Shape[1]*Y[1]+Shape[2]*Y[2]+Shape[3]*Y[3];
  M[2] = Shape[0]*Z[0]+Shape[1]*Z[1]+Shape[2]*Z[2]+Shape[3]*Z[3];
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::GetIsoParamMappingNormalJacobianProductTemp(Scalar* JNormal, Scalar* m, CoordSetT& cs)
{
  // Compute dM/dx & dM/dy
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // JN = dM/dx x dM/dy
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

#endif
