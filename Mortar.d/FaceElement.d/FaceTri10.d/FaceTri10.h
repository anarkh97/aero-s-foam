// ----------------------------------------------------------------
// HB - 05/05/05
// ----------------------------------------------------------------
#ifndef _FACETRI10_H_
#define _FACETRI10_H_

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

class FaceTri10: public FaceElement {
  private:
        int Nodes[10];
        static double RefCoords[10][2]; // coords of the nodes in the ref./parametric domain

  public:
        enum { NumberOfNodes=10 };

        // Constructors
        // ~~~~~~~~~~~~
        FaceTri10(int *);

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual methods
        void Renumber(std::map<int,int>& OldToNewNodeIds);

        // Get methods
        // ~~~~~~~~~~~
        // -> local methods
        int  nTri3Nodes();
        int  GetTri3Node(int);
        void GetTri3Nodes(int*, int* renumTable=0);
        void GetTri3Nodes(int*, std::map<int,int>& renumTable);

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
        template<typename Scalar>
          void   GetShapeFctVal(Scalar*, Scalar*);
        template<typename Scalar>
          void   GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar>
          void   Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m);
        template<typename Scalar>
          void   Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d2Shapex2y, Scalar *d2Shapexy2, Scalar *m);
        template<typename Scalar, typename CoordSetT>
          void   ComputedMdxAnddMdy(Scalar* dMdx, Scalar* dMdy, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          Scalar GetJacobian(Scalar*, CoordSetT&);
        template<typename Scalar, typename CoordSetT>
          Scalar GetShapeFctAndJacobian(Scalar* Shape, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   LocalToGlobalCoord(Scalar* M, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetIsoParamMappingNormalJacobianProduct(Scalar* JNormal, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          Scalar GetIsoParamMappingNormalAndJacobian(Scalar* Normal, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void Getd2JNormal(Scalar d2JNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs);

        // -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double*, double*, CoordSet&);
        void   GetShapeFctVal(double*, double*);
        double GetJacobian(double*, CoordSet&);
        double GetIsoParamMappingNormalAndJacobian(double*, double*, CoordSet&);
        void   GetIsoParamMappingNormalJacobianProduct(double*, double*, CoordSet&);

        // -> implementation of virtual fcts
        double* ViewRefCoords();
        void GetdShapeFct(double* dShapex, double* dShapey, double* m);
        void Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m);
        void Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m);
        void ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs);
        void Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs);
        void Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs);
        void GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs);
        void Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs);
        void ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs);
        void Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs);
        void ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs);

        // Miscelleaneous methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of virtual fcts
        //double ComputeArea(CoordSet&, const int ngp=2);

        // Mass matrix methods
        // ~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual fcts
        FullM ScalarMass(CoordSet&, double rho=1.0, int ngp=2);
        void  IntegrateShapeFcts(double*, CoordSet&, double rho=1.0, int ngp=2);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> local fcts
        void printNodes();

        // -> implementation of pure virtual fcts
        void print();

};

// -----------------------------------------------------------------------------------------------------
//                          MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
template<typename Scalar>
void
FaceTri10::GetShapeFctVal(Scalar *Shape, Scalar *m)
{
  // Computes the values of the shape functions at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P.
  // Outputs: Shape = [N₀(ξ,η), N₁(ξ,η), N₂(ξ,η), ⋯ , N₉(ξ,η)], the values of the shape functions at P.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s; 
  Scalar ninehalf= 9/2.;
  Scalar twonine = 2/9.;

  // following ACME Tri3 ordering ...
  Shape[0] = ninehalf*r*(r*(r-1)+twonine);
  Shape[1] = ninehalf*s*(s*(s-1)+twonine);
  Shape[2] = ninehalf*t*(t*(t-1)+twonine);
  Shape[3] = ninehalf*s*r*(3*r-1); 
  Shape[4] = ninehalf*s*r*(3*s-1);
  Shape[5] = ninehalf*t*s*(3*s-1); 
  Shape[6] = ninehalf*t*s*(3*t-1);
  Shape[7] = ninehalf*t*r*(3*t-1);
  Shape[8] = ninehalf*t*r*(3*r-1); 
  Shape[9] = 27*r*s*t;
}

template<typename Scalar>
void
FaceTri10::GetdShapeFct(Scalar *dShapex, Scalar *dShapey, Scalar *m)
{
  // Computes the values of the shape functions' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m       = [ξ, η], the local coordinates of P.
  // Outputs: dShapex = [∂N₀/∂ξ, ∂N₁/∂ξ, ∂N₂/∂ξ, ⋯ , ∂N₉/∂ξ]
  //          dShapey = [∂N₀/∂η, ∂N₁/∂η, ∂N₂/∂η, ⋯ , ∂N₉/∂η]
  //                    where N₀, N₁, N₂, ⋯ , and N₉ are the shape functions.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s;
  Scalar ninehalf= 9/2.;
  Scalar twonine = 2/9.;

  // !! idem ACME !!
  dShapex[0] =  ninehalf*(r*(3*r-2)+twonine);
  dShapex[1] =  0;
  dShapex[2] = -ninehalf*(t*(3*t-2)+twonine);
  dShapex[3] =  ninehalf*s*(6*r-1);
  dShapex[4] =  ninehalf*s*(3*s-1);
  dShapex[5] = -ninehalf*s*(3*s-1);
  dShapex[6] = -ninehalf*s*(6*t-1);
  dShapex[7] =  ninehalf*( t*(3*t-1) - r*(6*t-1));
  dShapex[8] =  ninehalf*(-r*(3*r-1) + t*(6*r-1));
  dShapex[9] =  27*s*(t-r);

  dShapey[0] =  0;
  dShapey[1] =  ninehalf*(s*(3*s-2)+twonine);
  dShapey[2] = -ninehalf*(t*(3*t-2)+twonine);
  dShapey[3] =  ninehalf*r*(3*r-1);
  dShapey[4] =  ninehalf*r*(6*s-1);
  dShapey[5] =  ninehalf*(-s*(3*s-1) + t*(6*s-1));
  dShapey[6] =  ninehalf*( t*(3*t-1) - s*(6*t-1));
  dShapey[7] = -ninehalf*r*(6*t-1);
  dShapey[8] = -ninehalf*r*(3*r-1);
  dShapey[9] =  27*r*(t-s);
}

template<typename Scalar>
void
FaceTri10::Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m)
{
  // Computes the values of shape functions' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m         = [ξ, η], the local coordinates of P.
  // Outputs: d2Shapex  = [∂²N₀/∂ξ², ∂²N₁/∂ξ², ∂²N₂/∂ξ², ⋯ , ∂²N₉/∂ξ²]
  //          d2Shapey  = [∂²N₀/∂η², ∂²N₁/∂η², ∂²N₂/∂η², ⋯ , ∂²N₉/∂η²]
  //          d2Shapexy = [∂²N₀/∂ξ∂η, ∂²N₁/∂ξ∂η, ∂²N₂/∂ξ∂η, ⋯ , ∂²N₉/∂ξ∂η]
  //                      where N₀, N₁, N₂, ⋯ , and N₉ are the shape functions.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s;
  Scalar ninehalf= 9/2.;
  Scalar twonine = 2/9.;

  d2Shapex[0] = ninehalf*(6*r-2);
  d2Shapex[1] = 0;
  d2Shapex[2] = ninehalf*(6*t-2);
  d2Shapex[3] = ninehalf*(6*s);
  d2Shapex[4] = 0;
  d2Shapex[5] = 0;
  d2Shapex[6] = ninehalf*(6*s);
  d2Shapex[7] = ninehalf*(-12*t+6*r+2);
  d2Shapex[8] = ninehalf*(-12*r+6*t+2);
  d2Shapex[9] = -54*s;

  d2Shapey[0] = 0;
  d2Shapey[1] = ninehalf*(6*s-2);
  d2Shapey[2] = ninehalf*(6*t-2);
  d2Shapey[3] = 0;
  d2Shapey[4] = ninehalf*(6*r);
  d2Shapey[5] = ninehalf*(-12*s+6*t+2);
  d2Shapey[6] = ninehalf*(-12*t+6*s+2);
  d2Shapey[7] = ninehalf*(6*r);
  d2Shapey[8] = 0;
  d2Shapey[9] = -54*r;

  d2Shapexy[0] = 0;
  d2Shapexy[1] = 0;
  d2Shapexy[2] = ninehalf*(6*t-2);
  d2Shapexy[3] = ninehalf*(6*r-1);
  d2Shapexy[4] = ninehalf*(6*s-1);
  d2Shapexy[5] = -ninehalf*(6*s-1);
  d2Shapexy[6] = -ninehalf*(-6*s+6*t-1);
  d2Shapexy[7] = ninehalf*(-6*t+6*r+1);
  d2Shapexy[8] = ninehalf*(-6*r+1);
  d2Shapexy[9] = 27*(-s+t-r);
}

template<typename Scalar>
void
FaceTri10::Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d3Shapex2y, Scalar *d3Shapexy2, Scalar *m)
{
  // Computes the values of the shape functions' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of P.
  // Outputs: d3Shapex   = [∂³N₀/∂ξ³, ∂³N₁/∂ξ³, ∂³N₂/∂ξ³, ⋯ , ∂³N₉/∂ξ³]
  //          d3Shapey   = [∂³N₀/∂η³, ∂³N₁/∂η³, ∂³N₂/∂η³, ⋯ , ∂³N₉/∂η³]
  //          d3Shapex2y = [∂³N₀/∂ξ²∂η, ∂³N₁/∂ξ²∂η, ∂³N₂/∂ξ²∂η, ⋯ , ∂³N₉/∂ξ²∂η]
  //          d3Shapexy2 = [∂³N₀/∂ξ∂η², ∂³N₁/∂ξ∂η², ∂³N₂/∂ξ∂η², ⋯ , ∂³N₉/∂ξ∂η²]
  //                       where N₀, N₁, N₂, ⋯ , and N₉ are the shape functions.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s;
  Scalar ninehalf= 9/2.;
  Scalar twonine = 2/9.;

  d3Shapex[0] = 27;
  d3Shapex[1] = 0;
  d3Shapex[2] = -27;
  d3Shapex[3] = 0;
  d3Shapex[4] = 0;
  d3Shapex[5] = 0;
  d3Shapex[6] = 0;
  d3Shapex[7] = 81;
  d3Shapex[8] = -81;
  d3Shapex[9] = 0;

  d3Shapey[0] = 0;
  d3Shapey[1] = 27;
  d3Shapey[2] = -27;
  d3Shapey[3] = 0;
  d3Shapey[4] = 0;
  d3Shapey[5] = -81;
  d3Shapey[6] = 81;
  d3Shapey[7] = 0;
  d3Shapey[8] = 0;
  d3Shapey[9] = 0;

  d3Shapex2y[0] = 0;
  d3Shapex2y[1] = 0;
  d3Shapex2y[2] = -27;
  d3Shapex2y[3] = 27;
  d3Shapex2y[4] = 0;
  d3Shapex2y[5] = 0;
  d3Shapex2y[6] = 27;
  d3Shapex2y[7] = 54;
  d3Shapex2y[8] = -27;
  d3Shapex2y[9] = -54;

  d3Shapexy2[0] = 0;
  d3Shapexy2[1] = 0;
  d3Shapexy2[2] = -27;
  d3Shapexy2[3] = 0;
  d3Shapexy2[4] = 27;
  d3Shapexy2[5] = -27;
  d3Shapexy2[6] = 54;
  d3Shapexy2[7] = 27;
  d3Shapexy2[8] = 0;
  d3Shapexy2[9] = -54;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::LocalToGlobalCoord(Scalar *M, Scalar *m, CoordSetT &cs)
{
  // Computes the global X-, Y- and Z-coordinates of a point P.
  //
  // Inputs:  m  = [ξ, η], the local coordinates of P, and
  //          cs = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: M  = [M₀(ξ,η,X,Y,Z), M₁(ξ,η,X,Y,Z), M₂(ξ,η,X,Y,Z)], the global X-, Y- and Z-coordinates of P.

  Scalar Shape[10];
  GetShapeFctVal(Shape,m);

  Scalar X[10], Y[10], Z[10];
  for(int i=0; i<10; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  M[0] = 0.0; M[1] = 0.0; M[2] = 0.0;
  for(int i=0; i<10; i+=5) {
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2] + Shape[i+3]*X[i+3] + Shape[i+4]*X[i+4];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2] + Shape[i+3]*Y[i+3] + Shape[i+4]*Y[i+4];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2] + Shape[i+3]*Z[i+3] + Shape[i+4]*Z[i+4];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dMdx  = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ]
  //          dMdy  = [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η]
  //                  where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[10], dShapey[10];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar X[10], Y[10], Z[10];
  for(int i=0; i<10; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<10; i+=5) {
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2] + dShapex[i+3]*X[i+3] + dShapex[i+4]*X[i+4];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2] + dShapex[i+3]*Y[i+3] + dShapex[i+4]*Y[i+4];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2] + dShapex[i+3]*Z[i+3] + dShapex[i+4]*Z[i+4];

    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2] + dShapey[i+3]*X[i+3] + dShapey[i+4]*X[i+4];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2] + dShapey[i+3]*Y[i+3] + dShapey[i+4]*Y[i+4];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2] + dShapey[i+3]*Z[i+3] + dShapey[i+4]*Z[i+4];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs: m        = [ξ, η], the local coordinates of P, and
  //         cs       = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2Mdx2  = [∂²M₀/∂ξ², ∂²M₁/∂ξ², ∂²M₂/∂ξ²]
  //          d2Mdy2  = [∂²M₀/∂η², ∂²M₁/∂η², ∂²M₂/∂η²]
  //          d2Mdxdy = [∂²M₀/∂ξ∂η, ∂²M₁/∂ξ∂η, ∂²M₂/∂ξ∂η]
  //                    where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Compute shape functions' 2nd derivatives w.r.t. the local coordinates  
  Scalar d2Shapex[10], d2Shapey[10], d2Shapexy[10];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar X[10], Y[10], Z[10];
  for(int i=0; i<10; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d2Mdx2[0] = d2Mdx2[1] = d2Mdx2[2] = 0.0;
  d2Mdy2[0] = d2Mdy2[1] = d2Mdy2[2] = 0.0;
  d2Mdxdy[0] = d2Mdxdy[1] = d2Mdxdy[2] = 0.0;
  for(int i=0; i<10; i+=5) {
    d2Mdx2[0] += d2Shapex[i]*X[i] + d2Shapex[i+1]*X[i+1] + d2Shapex[i+2]*X[i+2] + d2Shapex[i+3]*X[i+3] + d2Shapex[i+4]*X[i+4];
    d2Mdx2[1] += d2Shapex[i]*Y[i] + d2Shapex[i+1]*Y[i+1] + d2Shapex[i+2]*Y[i+2] + d2Shapex[i+3]*Y[i+3] + d2Shapex[i+4]*Y[i+4];
    d2Mdx2[2] += d2Shapex[i]*Z[i] + d2Shapex[i+1]*Z[i+1] + d2Shapex[i+2]*Z[i+2] + d2Shapex[i+3]*Z[i+3] + d2Shapex[i+4]*Z[i+4];

    d2Mdy2[0] += d2Shapey[i]*X[i] + d2Shapey[i+1]*X[i+1] + d2Shapey[i+2]*X[i+2] + d2Shapey[i+3]*X[i+3] + d2Shapey[i+4]*X[i+4];
    d2Mdy2[1] += d2Shapey[i]*Y[i] + d2Shapey[i+1]*Y[i+1] + d2Shapey[i+2]*Y[i+2] + d2Shapey[i+3]*Y[i+3] + d2Shapey[i+4]*Y[i+4];
    d2Mdy2[2] += d2Shapey[i]*Z[i] + d2Shapey[i+1]*Z[i+1] + d2Shapey[i+2]*Z[i+2] + d2Shapey[i+3]*Z[i+3] + d2Shapey[i+4]*Z[i+4];

    d2Mdxdy[0] += d2Shapexy[i]*X[i] + d2Shapexy[i+1]*X[i+1] + d2Shapexy[i+2]*X[i+2] + d2Shapexy[i+3]*X[i+3] + d2Shapexy[i+4]*X[i+4];
    d2Mdxdy[1] += d2Shapexy[i]*Y[i] + d2Shapexy[i+1]*Y[i+1] + d2Shapexy[i+2]*Y[i+2] + d2Shapexy[i+3]*Y[i+3] + d2Shapexy[i+4]*Y[i+4];
    d2Mdxdy[2] += d2Shapexy[i]*Z[i] + d2Shapexy[i+1]*Z[i+1] + d2Shapexy[i+2]*Z[i+2] + d2Shapexy[i+3]*Z[i+3] + d2Shapexy[i+4]*Z[i+4];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs: m         = [ξ, η], the local coordinates of P, and
  //         cs        = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d3Mdx3   = [∂³M₀/∂ξ³, ∂³M₁/∂ξ³, ∂³M₂/∂ξ³]
  //          d3Mdy3   = [∂³M₀/∂η³, ∂³M₁/∂η³, ∂³M₂/∂η³]
  //          d3Mdx2dy = [∂³M₀/∂ξ²∂η, ∂³M₁/∂ξ²∂η, ∂³M₂/∂ξ²∂η]
  //          d3Mdxdy2 = [∂³M₀/∂ξ∂η², ∂³M₁/∂ξ∂η², ∂³M₂/∂ξ∂η²]
  //                     where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Compute shape functions' 3rd derivatives w.r.t. the local coordinates  
  Scalar d3Shapex[10], d3Shapey[10], d3Shapex2y[10], d3Shapexy2[10];
  Getd3ShapeFct(d3Shapex, d3Shapey, d3Shapex2y, d3Shapexy2, m);

  // Compute ∂³M/∂ξ³, ∂³M/∂η³, ∂³M/∂ξ²∂η and ∂³M/∂ξ∂η²
  Scalar X[10], Y[10], Z[10];
  for(int i=0; i<10; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d3Mdx3[0] = d3Mdx3[1] = d3Mdx3[2] = 0.0;
  d3Mdy3[0] = d3Mdy3[1] = d3Mdy3[2] = 0.0;
  d3Mdx2dy[0] = d3Mdx2dy[1] = d3Mdx2dy[2] = 0.0;
  d3Mdxdy2[0] = d3Mdxdy2[1] = d3Mdxdy2[2] = 0.0;
  for(int i=0; i<10; i+=5) {
    d3Mdx3[0] += d3Shapex[i]*X[i] + d3Shapex[i+1]*X[i+1] + d3Shapex[i+2]*X[i+2] + d3Shapex[i+3]*X[i+3] + d3Shapex[i+4]*X[i+4];
    d3Mdx3[1] += d3Shapex[i]*Y[i] + d3Shapex[i+1]*Y[i+1] + d3Shapex[i+2]*Y[i+2] + d3Shapex[i+3]*Y[i+3] + d3Shapex[i+4]*Y[i+4];
    d3Mdx3[2] += d3Shapex[i]*Z[i] + d3Shapex[i+1]*Z[i+1] + d3Shapex[i+2]*Z[i+2] + d3Shapex[i+3]*Z[i+3] + d3Shapex[i+4]*Z[i+4];

    d3Mdy3[0] += d3Shapey[i]*X[i] + d3Shapey[i+1]*X[i+1] + d3Shapey[i+2]*X[i+2] + d3Shapey[i+3]*X[i+3] + d3Shapey[i+4]*X[i+4];
    d3Mdy3[1] += d3Shapey[i]*Y[i] + d3Shapey[i+1]*Y[i+1] + d3Shapey[i+2]*Y[i+2] + d3Shapey[i+3]*Y[i+3] + d3Shapey[i+4]*Y[i+4];
    d3Mdy3[2] += d3Shapey[i]*Z[i] + d3Shapey[i+1]*Z[i+1] + d3Shapey[i+2]*Z[i+2] + d3Shapey[i+3]*Z[i+3] + d3Shapey[i+4]*Z[i+4];

    d3Mdx2dy[0] += d3Shapex2y[i]*X[i] + d3Shapex2y[i+1]*X[i+1] + d3Shapex2y[i+2]*X[i+2] + d3Shapex2y[i+3]*X[i+3] + d3Shapex2y[i+4]*X[i+4];
    d3Mdx2dy[1] += d3Shapex2y[i]*Y[i] + d3Shapex2y[i+1]*Y[i+1] + d3Shapex2y[i+2]*Y[i+2] + d3Shapex2y[i+3]*Y[i+3] + d3Shapex2y[i+4]*Y[i+4];
    d3Mdx2dy[2] += d3Shapex2y[i]*Z[i] + d3Shapex2y[i+1]*Z[i+1] + d3Shapex2y[i+2]*Z[i+2] + d3Shapex2y[i+3]*Z[i+3] + d3Shapex2y[i+4]*Z[i+4];

    d3Mdxdy2[0] += d3Shapexy2[i]*X[i] + d3Shapexy2[i+1]*X[i+1] + d3Shapexy2[i+2]*X[i+2] + d3Shapexy2[i+3]*X[i+3] + d3Shapexy2[i+4]*X[i+4];
    d3Mdxdy2[1] += d3Shapexy2[i]*Y[i] + d3Shapexy2[i+1]*Y[i+1] + d3Shapexy2[i+2]*Y[i+2] + d3Shapexy2[i+3]*Y[i+3] + d3Shapexy2[i+4]*Y[i+4];
    d3Mdxdy2[2] += d3Shapexy2[i]*Z[i] + d3Shapexy2[i+1]*Z[i+1] + d3Shapexy2[i+2]*Z[i+2] + d3Shapexy2[i+3]*Z[i+3] + d3Shapexy2[i+4]*Z[i+4];
  }
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri10::GetJacobian(Scalar *m, CoordSetT &cs)
{
  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // JNormal = ∂M/∂ξ x ∂M/∂η
  Scalar JNormal[3];
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  using std::sqrt;
  return(sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]));
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri10::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  GetShapeFctVal(Shape, m);
  return(GetJacobian(m, cs));
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::GetIsoParamMappingNormalJacobianProduct(Scalar *JNormal, Scalar *m, CoordSetT &cs)
{
  // Computes the vector normal to the surface at a point P whose magnitude is the Jacobian determinant
  // of the local-to-global mapping.

  // Inputs:  m       = [ξ, η], the local coordinates of P, and
  //          cs      = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: JNormal = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ] x [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η], the vector normal
  //                    to the surface at point P whose magnitude is the Jacobian determinant of
  //                    the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                    global X-, Y- and Z-coordinates of P.

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute the normal vector, n = ∂M/∂ξ x ∂M/∂η
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the partial derivatives w.r.t the local coordinates of the components of the vector
  // normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of some point P, and
  //          cs         = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormaldx = [∂n₀/∂ξ, ∂n₁/∂ξ, ∂n₂/∂ξ]
  //          dJNormaldy = [∂n₀/∂η, ∂n₁/∂η, ∂n₂/∂η]
  //                       where n₀, n₁, and n₂ are the global X-, Y- and Z-components of the vector
  //                       normal to the surface at the point P whose magnitude is the Jacobian determinant
  //                       of the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                       global X-, Y- and Z-coordinates of P.

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂n/∂ξ and ∂n/∂η
  dJNormaldx[0] = d2Mdx2[1]*dMdy[2] + dMdx[1]*d2Mdxdy[2] - (d2Mdx2[2]*dMdy[1] + dMdx[2]*d2Mdxdy[1]);
  dJNormaldx[1] = d2Mdx2[2]*dMdy[0] + dMdx[2]*d2Mdxdy[0] - (d2Mdx2[0]*dMdy[2] + dMdx[0]*d2Mdxdy[2]);
  dJNormaldx[2] = d2Mdx2[0]*dMdy[1] + dMdx[0]*d2Mdxdy[1] - (d2Mdx2[1]*dMdy[0] + dMdx[1]*d2Mdxdy[0]);

  dJNormaldy[0] = d2Mdxdy[1]*dMdy[2] + dMdx[1]*d2Mdy2[2] - (d2Mdxdy[2]*dMdy[1] + dMdx[2]*d2Mdy2[1]);
  dJNormaldy[1] = d2Mdxdy[2]*dMdy[0] + dMdx[2]*d2Mdy2[0] - (d2Mdxdy[0]*dMdy[2] + dMdx[0]*d2Mdy2[2]);
  dJNormaldy[2] = d2Mdxdy[0]*dMdy[1] + dMdx[0]*d2Mdy2[1] - (d2Mdxdy[1]*dMdy[0] + dMdx[1]*d2Mdy2[0]);
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the second partial derivatives w.r.t the local coordinates of the components of the vector
  // normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m             = [ξ, η], the local coordinates of some point P, and
  //          cs            = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2JNormaldx2  = [∂²n₀/∂ξ², ∂²n₁/∂ξ², ∂²n₂/∂ξ²]
  //          d2JNormaldy2  = [∂²n₀/∂η², ∂²n₁/∂η², ∂²n₂/∂η²]
  //          d2JNormaldxdy = [∂²n₀/∂ξ∂η, ∂²n₁/∂ξ∂η, ∂²n₂/∂ξ∂η]
  //                          where n₀, n₁, and n₂ are the global X-, Y- and Z-components of the vector
  //                          normal to the surface at the point P whose magnitude is the Jacobian determinant
  //                          of the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                          global X-, Y- and Z-coordinates of P.

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂³M/∂ξ³, ∂³M/∂η³, ∂³M/∂ξ²∂η and ∂³M/∂ξ∂η²
  Scalar d3Mdx3[3], d3Mdy3[3], d3Mdx2dy[3], d3Mdxdy2[3];
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);

  // Compute ∂²n/∂ξ², ∂²n/∂η² and ∂²n/∂ξ∂η
  d2JNormaldx2[0] = d2Mdx2[1]*d2Mdxdy[2] + d2Mdx2[1]*d2Mdxdy[2] - (d2Mdx2[2]*d2Mdxdy[1] + d2Mdx2[2]*d2Mdxdy[1])
                  + d3Mdx3[1]*dMdy[2] + dMdx[1]*d3Mdx2dy[2] - (d3Mdx3[2]*dMdy[1] + dMdx[2]*d3Mdx2dy[1]);
  d2JNormaldx2[1] = d2Mdx2[2]*d2Mdxdy[0] + d2Mdx2[2]*d2Mdxdy[0] - (d2Mdx2[0]*d2Mdxdy[2] + d2Mdx2[0]*d2Mdxdy[2])
                  + d3Mdx3[2]*dMdy[0] + dMdx[2]*d3Mdx2dy[0] - (d3Mdx3[0]*dMdy[2] + dMdx[0]*d3Mdx2dy[2]);
  d2JNormaldx2[2] = d2Mdx2[0]*d2Mdxdy[1] + d2Mdx2[0]*d2Mdxdy[1] - (d2Mdx2[1]*d2Mdxdy[0] + d2Mdx2[1]*d2Mdxdy[0])
                  + d3Mdx3[0]*dMdy[1] + dMdx[0]*d3Mdx2dy[1] - (d3Mdx3[1]*dMdy[0] + dMdx[1]*d3Mdx2dy[0]);

  d2JNormaldy2[0] = d2Mdxdy[1]*d2Mdy2[2] + d2Mdxdy[1]*d2Mdy2[2] - (d2Mdxdy[2]*d2Mdy2[1] + d2Mdxdy[2]*d2Mdy2[1])
                  + d3Mdxdy2[1]*dMdy[2] + dMdx[1]*d3Mdy3[2] - (d3Mdxdy2[2]*dMdy[1] + dMdx[2]*d3Mdy3[1]);
  d2JNormaldy2[1] = d2Mdxdy[2]*d2Mdy2[0] + d2Mdxdy[2]*d2Mdy2[0] - (d2Mdxdy[0]*d2Mdy2[2] + d2Mdxdy[0]*d2Mdy2[2])
                  + d3Mdxdy2[2]*dMdy[0] + dMdx[2]*d3Mdy3[0] - (d3Mdxdy2[0]*dMdy[2] + dMdx[0]*d3Mdy3[2]);
  d2JNormaldy2[2] = d2Mdxdy[0]*d2Mdy2[1] + d2Mdxdy[0]*d2Mdy2[1] - (d2Mdxdy[1]*d2Mdy2[0] + d2Mdxdy[1]*d2Mdy2[0])
                  + d3Mdxdy2[0]*dMdy[1] + dMdx[0]*d3Mdy3[1] - (d3Mdxdy2[1]*dMdy[0] + dMdx[1]*d3Mdy3[0]);

  d2JNormaldxdy[0] = d2Mdx2[1]*d2Mdy2[2] + d2Mdxdy[1]*d2Mdxdy[2] - (d2Mdx2[2]*d2Mdy2[1] + d2Mdxdy[2]*d2Mdxdy[1])
                   + d3Mdx2dy[1]*dMdy[2] + dMdx[1]*d3Mdxdy2[2] - (d3Mdx2dy[2]*dMdy[1] + dMdx[2]*d3Mdxdy2[1]);
  d2JNormaldxdy[1] = d2Mdx2[2]*d2Mdy2[0] + d2Mdxdy[2]*d2Mdxdy[0] - (d2Mdx2[0]*d2Mdy2[2] + d2Mdxdy[0]*d2Mdxdy[2])
                   + d3Mdx2dy[2]*dMdy[0] + dMdx[2]*d3Mdxdy2[0] - (d3Mdx2dy[0]*dMdy[2] + dMdx[0]*d3Mdxdy2[2]);
  d2JNormaldxdy[2] = d2Mdx2[0]*d2Mdy2[1] + d2Mdxdy[0]*d2Mdxdy[1] - (d2Mdx2[1]*d2Mdy2[0] + d2Mdxdy[1]*d2Mdxdy[0])
                   + d3Mdx2dy[0]*dMdy[1] + dMdx[0]*d3Mdxdy2[1] - (d3Mdx2dy[1]*dMdy[0] + dMdx[1]*d3Mdxdy2[0]);
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri10::GetIsoParamMappingNormalAndJacobian(Scalar *Normal, Scalar *m, CoordSetT &cs)
{
  // Compute dM/dx and dM/dy
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // JN = dM/dx x dM/dy
  Normal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  Normal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  Normal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  using std::sqrt;
  Scalar NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(NormN != 0.0) {
    Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
  }
  return(NormN);
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.

  // Inputs:  m           = [ξ, η], the local coordinates of P, and
  //          cs          = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormal[0] = [∂n₀/∂X₀, ∂n₀/∂Y₀, ∂n₀/∂Z₀, ∂n₀/∂X₁, ∂n₀/∂Y₁, ∂n₀/∂Z₁, ⋯ , ∂n₀/∂X₉, ∂n₀/∂Y₉, ∂n₀/∂Z₉]
  //          dJNormal[1] = [∂n₁/∂X₀, ∂n₁/∂Y₀, ∂n₁/∂Z₀, ∂n₁/∂X₁, ∂n₁/∂Y₁, ∂n₁/∂Z₁, ⋯ , ∂n₁/∂X₉, ∂n₁/∂Y₉, ∂n₁/∂Z₉]
  //          dJNormal[2] = [∂n₂/∂X₀, ∂n₂/∂Y₀, ∂n₂/∂Z₀, ∂n₂/∂X₁, ∂n₂/∂Y₁, ∂n₂/∂Z₁, ⋯ , ∂n₂/∂X₉, ∂n₂/∂Y₉, ∂n₂/∂Z₉]
  //                        where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                        at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                        where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                        and X = [X₀, X₁, ⋯ , X₉], Y = [Y₀, Y₁, ⋯ , Y₉] and  Z = [Z₀, Z₁, ⋯ , Z₉] are
  //                        the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[10], dShapey[10];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute dJNormal
  for(int i = 0; i < 10; ++i) {
    dJNormal[3*i  ][0] = 0;
    dJNormal[3*i  ][1] = dMdx[2]*dShapey[i] - dShapex[i]*dMdy[2];
    dJNormal[3*i  ][2] = dShapex[i]*dMdy[1] - dMdx[1]*dShapey[i];

    dJNormal[3*i+1][0] = dShapex[i]*dMdy[2] - dMdx[2]*dShapey[i];
    dJNormal[3*i+1][1] = 0;
    dJNormal[3*i+1][2] = dMdx[0]*dShapey[i] - dShapex[i]*dMdy[0];

    dJNormal[3*i+2][0] = dMdx[1]*dShapey[i] - dShapex[i]*dMdy[1];
    dJNormal[3*i+2][1] = dShapex[i]*dMdy[0] - dMdx[0]*dShapey[i];
    dJNormal[3*i+2][2] = 0;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the mixed second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates and the local coordinate
  // of the components of the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m                = [ξ, η], the local coordinates of P, and
  //          cs               = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: ddJNormaldx[][0] = [∂²n₀/∂ξ∂X₀, ∂²n₀/∂ξ∂Y₀, ∂²n₀/∂ξ∂Z₀, ∂²n₀/∂ξ∂X₁, ∂²n₀/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₀/∂ξ∂X₉, ∂²n₀/∂ξ∂Y₉, ∂²n₀/∂ξ∂Z₉]
  //          ddJNormaldx[][1] = [∂²n₁/∂ξ∂X₀, ∂²n₁/∂ξ∂Y₀, ∂²n₁/∂ξ∂Z₀, ∂²n₁/∂ξ∂X₁, ∂²n₁/∂ξ∂Y₁, ∂²n₁/∂ξ∂Z₁, ⋯ , ∂²n₁/∂ξ∂X₉, ∂²n₁/∂ξ∂Y₉, ∂²n₁/∂ξ∂Z₉]
  //          ddJNormaldx[][2] = [∂²n₂/∂ξ∂X₀, ∂²n₂/∂ξ∂Y₀, ∂²n₂/∂ξ∂Z₀, ∂²n₂/∂ξ∂X₁, ∂²n₂/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₂/∂ξ∂X₉, ∂²n₂/∂ξ∂Y₉, ∂²n₂/∂ξ∂Z₉]
  //          ddJNormaldy[][0] = [∂²n₀/∂η∂X₀, ∂²n₀/∂η∂Y₀, ∂²n₀/∂η∂Z₀, ∂²n₀/∂η∂X₁, ∂²n₀/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₀/∂η∂X₉, ∂²n₀/∂η∂Y₉, ∂²n₀/∂η∂Z₉]
  //          ddJNormaldy[][1] = [∂²n₁/∂η∂X₀, ∂²n₁/∂η∂Y₀, ∂²n₁/∂η∂Z₀, ∂²n₁/∂η∂X₁, ∂²n₁/∂η∂Y₁, ∂²n₁/∂η∂Z₁, ⋯ , ∂²n₁/∂η∂X₉, ∂²n₁/∂η∂Y₉, ∂²n₁/∂η∂Z₉]
  //          ddJNormaldy[][2] = [∂²n₂/∂η∂X₀, ∂²n₂/∂η∂Y₀, ∂²n₂/∂η∂Z₀, ∂²n₂/∂η∂X₁, ∂²n₂/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₂/∂η∂X₉, ∂²n₂/∂η∂Y₉, ∂²n₂/∂η∂Z₉]
  //                             where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                             at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                             where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                             and X = [X₀, X₁, ⋯ , X₉], Y = [Y₀, Y₁, ⋯ , Y₉] and  Z = [Z₀, Z₁, ⋯ , Z₉] are
  //                             the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[10], dShapey[10];
  GetdShapeFct(dShapex, dShapey, m);
  
  // Compute shape functions' 2nd derivatives w.r.t. the local coordinates
  Scalar d2Shapex[10], d2Shapey[10], d2Shapexy[10];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ddJNormaldx and ddJNormaldy
  for(int i = 0; i < 10; ++i) {
    ddJNormaldx[3*i  ][0] = 0;
    ddJNormaldx[3*i  ][1] = (d2Mdx2[2]*dShapey[i] + dMdx[2]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[2] + dShapex[i]*d2Mdxdy[2]);
    ddJNormaldx[3*i  ][2] = (d2Shapex[i]*dMdy[1] + dShapex[i]*d2Mdxdy[1]) - (d2Mdx2[1]*dShapey[i] + dMdx[1]*d2Shapexy[i]);
    ddJNormaldx[3*i+1][0] = (d2Shapex[i]*dMdy[2] + dShapex[i]*d2Mdxdy[2]) - (d2Mdx2[2]*dShapey[i] + dMdx[2]*d2Shapexy[i]);
    ddJNormaldx[3*i+1][1] = 0;
    ddJNormaldx[3*i+1][2] = (d2Mdx2[0]*dShapey[i] + dMdx[0]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[0] + dShapex[i]*d2Mdxdy[0]);
    ddJNormaldx[3*i+2][0] = (d2Mdx2[1]*dShapey[i] + dMdx[1]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[1] + dShapex[i]*d2Mdxdy[1]);
    ddJNormaldx[3*i+2][1] = (d2Shapex[i]*dMdy[0] + dShapex[i]*d2Mdxdy[0]) - (d2Mdx2[0]*dShapey[i] + dMdx[0]*d2Shapexy[i]);
    ddJNormaldx[3*i+2][2] = 0;
  }

  for(int i = 0; i < 10; ++i) {
    ddJNormaldy[3*i  ][0] = 0;
    ddJNormaldy[3*i  ][1] = (d2Mdxdy[2]*dShapey[i] + dMdx[2]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[2] + dShapex[i]*d2Mdy2[2]);
    ddJNormaldy[3*i  ][2] = (d2Shapexy[i]*dMdy[1] + dShapex[i]*d2Mdy2[1]) - (d2Mdxdy[1]*dShapey[i] + dMdx[1]*d2Shapey[i]);
    ddJNormaldy[3*i+1][0] = (d2Shapexy[i]*dMdy[2] + dShapex[i]*d2Mdy2[2]) - (d2Mdxdy[2]*dShapey[i] + dMdx[2]*d2Shapey[i]);
    ddJNormaldy[3*i+1][1] = 0;
    ddJNormaldy[3*i+1][2] = (d2Mdxdy[0]*dShapey[i] + dMdx[0]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[0] + dShapex[i]*d2Mdy2[0]);
    ddJNormaldy[3*i+2][0] = (d2Mdxdy[1]*dShapey[i] + dMdx[1]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[1] + dShapex[i]*d2Mdy2[1]);
    ddJNormaldy[3*i+2][1] = (d2Shapexy[i]*dMdy[0] + dShapex[i]*d2Mdy2[0]) - (d2Mdxdy[0]*dShapey[i] + dMdx[0]*d2Shapey[i]);
    ddJNormaldy[3*i+2][2] = 0;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri10::Getd2JNormal(Scalar H[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m    = [ξ, η], the local coordinates of P, and
  //          cs   = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: H[0] = [∂²n₀/∂X₀²,   ∂²n₀/∂X₀∂Y₀, ∂²n₀/∂X₀∂Z₀, ∂²n₀/∂X₀∂X₁, ∂²n₀/∂X₀∂Y₁, ∂²n₀/∂X₀∂Z₁, ⋯ , ∂²n₀/∂X₀∂X₉, ∂²n₀/∂X₀∂Y₉, ∂²n₀/∂X₀∂Z₉,
  //                  ∂²n₀/∂Y₀∂X₀, ∂²n₀/∂Y₀²,   ∂²n₀/∂Y₀∂Z₀, ∂²n₀/∂Y₀∂X₁, ∂²n₀/∂Y₀∂Y₁, ∂²n₀/∂Y₀∂Z₁, ⋯ , ∂²n₀/∂Y₀∂X₉, ∂²n₀/∂Y₀∂Y₉, ∂²n₀/∂Y₀∂Z₉,
  //                                                                         ⋮ 
  //                  ∂²n₀/∂Z₉∂X₀, ∂²n₀/∂Z₉∂Y₀, ∂²n₀/∂Z₉∂Z₀, ∂²n₀/∂Z₉∂X₁, ∂²n₀/∂Z₉∂Y₁, ∂²n₀/∂Z₉∂Z₁, ⋯ , ∂²n₀/∂Z₉∂X₉, ∂²n₀/∂Z₉∂Y₉, ∂²n₀/∂Z₉²]
  //          H[1] = [∂²n₁/∂X₀²,   ∂²n₁/∂X₀∂Y₀, ∂²n₁/∂X₀∂Z₀, ∂²n₁/∂X₀∂X₁, ∂²n₁/∂X₀∂Y₁, ∂²n₁/∂X₀∂Z₁, ⋯ , ∂²n₁/∂X₀∂X₉, ∂²n₁/∂X₀∂Y₉, ∂²n₁/∂X₀∂Z₉,
  //                  ∂²n₁/∂Y₀∂X₀, ∂²n₁/∂Y₀²,   ∂²n₁/∂Y₀∂Z₀, ∂²n₁/∂Y₀∂X₁, ∂²n₁/∂Y₀∂Y₁, ∂²n₁/∂Y₀∂Z₁, ⋯ , ∂²n₁/∂Y₀∂X₉, ∂²n₁/∂Y₀∂Y₉, ∂²n₁/∂Y₀∂Z₉,
  //                                                                         ⋮ 
  //                  ∂²n₁/∂Z₉∂X₀, ∂²n₁/∂Z₉∂Y₀, ∂²n₁/∂Z₉∂Z₀, ∂²n₁/∂Z₉∂X₁, ∂²n₁/∂Z₉∂Y₁, ∂²n₁/∂Z₉∂Z₁, ⋯ , ∂²n₁/∂Z₉∂X₉, ∂²n₁/∂Z₉∂Y₉, ∂²n₁/∂Z₉²]
  //          H[2] = [∂²n₂/∂X₀²,   ∂²n₂/∂X₀∂Y₀, ∂²n₂/∂X₀∂Z₀, ∂²n₂/∂X₀∂X₁, ∂²n₂/∂X₀∂Y₁, ∂²n₂/∂X₀∂Z₁, ⋯ , ∂²n₂/∂X₀∂X₉, ∂²n₂/∂X₀∂Y₉, ∂²n₂/∂X₀∂Z₉,
  //                  ∂²n₂/∂Y₀∂X₀, ∂²n₂/∂Y₀²,   ∂²n₂/∂Y₀∂Z₀, ∂²n₂/∂Y₀∂X₁, ∂²n₂/∂Y₀∂Y₁, ∂²n₂/∂Y₀∂Z₁, ⋯ , ∂²n₂/∂Y₀∂X₉, ∂²n₂/∂Y₀∂Y₉, ∂²n₂/∂Y₀∂Z₉,
  //                                                                         ⋮ 
  //                  ∂²n₂/∂Z₉∂X₀, ∂²n₂/∂Z₉∂Y₀, ∂²n₂/∂Z₉∂Z₀, ∂²n₂/∂Z₉∂X₁, ∂²n₂/∂Z₉∂Y₁, ∂²n₂/∂Z₉∂Z₁, ⋯ , ∂²n₂/∂Z₉∂X₉, ∂²n₂/∂Z₉∂Y₉, ∂²n₂/∂Z₉²]
  //                 where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                 at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                 where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                 and X = [X₀, X₁, ⋯ , X₉], Y = [Y₀, Y₁, ⋯ , Y₉] and  Z = [Z₀, Z₁, ⋯ , Z₉] are
  //                 the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' first partial derivatives w.r.t. the local coordinates
  Scalar dShapex[10], dShapey[10];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute d2JNormal
  for(int i = 0; i < 10; ++ i) {
    for(int j = 0; j < 10; ++j) {
      H[30*(3*j  )+3*i  ][0] = 0;
      H[30*(3*j  )+3*i  ][1] = 0;
      H[30*(3*j  )+3*i  ][2] = 0;

      H[30*(3*j  )+3*i+1][0] = 0;
      H[30*(3*j  )+3*i+1][1] = 0;
      H[30*(3*j  )+3*i+1][2] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];

      H[30*(3*j  )+3*i+2][0] = 0;
      H[30*(3*j  )+3*i+2][1] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[30*(3*j  )+3*i+2][2] = 0;

      H[30*(3*j+1)+3*i  ][0] = 0;
      H[30*(3*j+1)+3*i  ][1] = 0;
      H[30*(3*j+1)+3*i  ][2] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];

      H[30*(3*j+1)+3*i+1][0] = 0;
      H[30*(3*j+1)+3*i+1][1] = 0;
      H[30*(3*j+1)+3*i+1][2] = 0;

      H[30*(3*j+1)+3*i+2][0] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[30*(3*j+1)+3*i+2][1] = 0;
      H[30*(3*j+1)+3*i+2][2] = 0;

      H[30*(3*j+2)+3*i  ][0] = 0;
      H[30*(3*j+2)+3*i  ][1] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[30*(3*j+2)+3*i  ][2] = 0;

      H[30*(3*j+2)+3*i+1][0] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[30*(3*j+2)+3*i+1][1] = 0;
      H[30*(3*j+2)+3*i+1][2] = 0;

      H[30*(3*j+2)+3*i+2][0] = 0;
      H[30*(3*j+2)+3*i+2][1] = 0;
      H[30*(3*j+2)+3*i+2][2] = 0;
    }
  }
}

#endif
