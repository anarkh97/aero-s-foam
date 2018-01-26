// ----------------------------------------------------------------
// HB - 08/26/03
// ----------------------------------------------------------------
#ifndef _FACEQUAD8_H_
#define _FACEQUAD8_H_

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

class FaceQuad8: public FaceElement {
  private:
        int Nodes[8];
        static double RefCoords[8][2]; // coords of the nodes in the ref./parametric domain

  public:
        enum { NumberOfNodes=8 };

        // Constructors
        // ~~~~~~~~~~~~
        FaceQuad8(int *);
        FaceElement* clone();

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
        int  nNodes() const override;
        void GetNodes(int*, int* renumTable) const override;
        void GetNodes(int*, std::map<int,int>& renumTable) const override;
        int  GetNode(int) const override;
        int  GetNodeIndex(int) const override;

        int GetFaceElemType() override;
#ifdef USE_ACME
        ContactSearch::ContactFace_Type GetACMEFaceElemType();
#else
        int GetACMEFaceElemType();
#endif
        // -> pure virtual method for dealing with quadratic face element
        //    (see FaceElement.h for more details)
        int  nVertices() override;
        int  GetVertex(int) override;
        void GetVertices(int*, int* renumTable) override;
        void GetVertices(int*, std::map<int,int>& renumTable) override;

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
        template<typename Scalar, typename CoordSetT>
          void GetUnitNormal(Scalar UnitNormal[3], Scalar* m, CoordSetT& cs);

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
        void GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs);

        // Miscelleaneous methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of virtual fcts
        //double ComputeArea(CoordSet&, const int ngp=3);

        // Mass matrix methods
        // ~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual fcts
        FullM ScalarMass(CoordSet&, double rho=1.0, int ngp=3);
        void  IntegrateShapeFcts(double*, CoordSet&, double rho=1.0, int ngp=3);

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
FaceQuad8::GetShapeFctVal(Scalar *Shape, Scalar *m)
{
  // Computes the values of the shape functions at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P.
  // Outputs: Shape = [N₀(ξ,η), N₁(ξ,η), N₂(ξ,η), ⋯ , N₇(ξ,η)], the values of the shape functions at P.

  Scalar& x = m[0];
  Scalar& y = m[1];
  Scalar onequart = 1/4.;
  Scalar onehalf  = 1/2.;

  Scalar xm = 1-x;
  Scalar xp = 1+x;
  Scalar ym = 1-y;
  Scalar yp = 1+y;
  Scalar xxm= 1-x*x;
  Scalar yym= 1-y*y;

  Shape[0] = onequart*xm*ym*(-1-x-y);
  Shape[1] = onequart*xp*ym*(-1+x-y);
  Shape[2] = onequart*xp*yp*(-1+x+y);
  Shape[3] = onequart*xm*yp*(-1-x+y);
  Shape[4] = onehalf*xxm*ym;
  Shape[5] = onehalf*yym*xp;
  Shape[6] = onehalf*xxm*yp;
  Shape[7] = onehalf*yym*xm;
}

template<typename Scalar>
void
FaceQuad8::GetdShapeFct(Scalar *dShapex, Scalar *dShapey, Scalar *m)
{
  // Computes the values of the shape functions' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m       = [ξ, η], the local coordinates of P.
  // Outputs: dShapex = [∂N₀/∂ξ, ∂N₁/∂ξ, ∂N₂/∂ξ, ⋯ , ∂N₇/∂ξ]
  //          dShapey = [∂N₀/∂η, ∂N₁/∂η, ∂N₂/∂η, ⋯ , ∂N₇/∂η]
  //                    where N₀, N₁, N₂, ⋯ , and N₇ are the shape functions.

  Scalar& x = m[0];
  Scalar& y = m[1];
  Scalar onequart = 1/4.;
  Scalar onehalf  = 1/2.;

  Scalar xm = 1-x;
  Scalar xp = 1+x;
  Scalar ym = 1-y;
  Scalar yp = 1+y;
  Scalar xxm= 1-x*x;
  Scalar yym= 1-y*y;

  dShapex[0] = onequart*ym*(2*x+y);
  dShapex[1] = onequart*ym*(2*x-y);
  dShapex[2] = onequart*yp*(2*x+y);
  dShapex[3] = onequart*yp*(2*x-y);
  dShapex[4] = -x*ym;
  dShapex[5] = onehalf*yym;
  dShapex[6] = -x*yp;
  dShapex[7] = -onehalf*yym;

  dShapey[0] = onequart*xm*( x+2*y);
  dShapey[1] = onequart*xp*(-x+2*y);
  dShapey[2] = onequart*xp*( x+2*y);
  dShapey[3] = onequart*xm*(-x+2*y);
  dShapey[4] = -onehalf*xxm;
  dShapey[5] = -y*xp;
  dShapey[6] = onehalf*xxm;
  dShapey[7] = -y*xm;
}

template<typename Scalar>
void
FaceQuad8::Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m)
{
  // Computes the values of the shape functions' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m         = [ξ, η], the local coordinates of P.
  // Outputs: d2Shapex  = [∂²N₀/∂ξ², ∂²N₁/∂ξ², ∂²N₂/∂ξ², ⋯ , ∂²N₇/∂ξ²]
  //          d2Shapey  = [∂²N₀/∂η², ∂²N₁/∂η², ∂²N₂/∂η², ⋯ , ∂²N₇/∂η²]
  //          d2Shapexy = [∂²N₀/∂ξ∂η, ∂²N₁/∂ξ∂η, ∂²N₂/∂ξ∂η, ⋯ , ∂²N₇/∂ξ∂η]
  //                      where N₀, N₁, N₂, ⋯ , and N₇ are the shape functions.

  Scalar& x = m[0];
  Scalar& y = m[1];
  Scalar onequart = 1/4.;
  Scalar onehalf  = 1/2.;

  Scalar xm = 1-x;
  Scalar xp = 1+x;
  Scalar ym = 1-y;
  Scalar yp = 1+y;

  d2Shapex[0] = onehalf*ym;
  d2Shapex[1] = onehalf*ym;
  d2Shapex[2] = onehalf*yp;
  d2Shapex[3] = onehalf*yp;
  d2Shapex[4] = -ym;
  d2Shapex[5] = 0;
  d2Shapex[6] = -yp;
  d2Shapex[7] = 0;

  d2Shapey[0] = onehalf*xm;
  d2Shapey[1] = onehalf*xp;
  d2Shapey[2] = onehalf*xp;
  d2Shapey[3] = onehalf*xm;
  d2Shapey[4] = 0;
  d2Shapey[5] = -xp;
  d2Shapey[6] = 0;
  d2Shapey[7] = -xm;

  d2Shapexy[0] = onequart*( 1-2*x-2*y);
  d2Shapexy[1] = onequart*(-1-2*x+2*y);
  d2Shapexy[2] = onequart*( 1+2*x+2*y);
  d2Shapexy[3] = onequart*(-1+2*x-2*y);
  d2Shapexy[4] = x;
  d2Shapexy[5] = -y;
  d2Shapexy[6] = -x;
  d2Shapexy[7] = y;
}

template<typename Scalar>
void
FaceQuad8::Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d3Shapex2y, Scalar *d3Shapexy2, Scalar *m)
{
  // Computes the values of the shape functions' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of P.
  // Outputs: d3Shapex   = [∂³N₀/∂ξ³, ∂³N₁/∂ξ³, ∂³N₂/∂ξ³, ⋯ , ∂³N₇/∂ξ³]
  //          d3Shapey   = [∂³N₀/∂η³, ∂³N₁/∂η³, ∂³N₂/∂η³, ⋯ , ∂³N₇/∂η³]
  //          d3Shapex2y = [∂³N₀/∂ξ²∂η, ∂³N₁/∂ξ²∂η, ∂³N₂/∂ξ²∂η, ⋯ , ∂³N₇/∂ξ²∂η]
  //          d3Shapexy2 = [∂³N₀/∂ξ∂η², ∂³N₁/∂ξ∂η², ∂³N₂/∂ξ∂η², ⋯ , ∂³N₇/∂ξ∂η²]
  //                       where N₀, N₁, N₂, ⋯ , and N₇ are the shape functions.

  Scalar onehalf  = 1/2.;

  d3Shapex[0] = 0;
  d3Shapex[1] = 0;
  d3Shapex[2] = 0;
  d3Shapex[3] = 0;
  d3Shapex[4] = 0;
  d3Shapex[5] = 0;
  d3Shapex[6] = 0;
  d3Shapex[7] = 0;

  d3Shapey[0] = 0;
  d3Shapey[1] = 0;
  d3Shapey[2] = 0;
  d3Shapey[3] = 0;
  d3Shapey[4] = 0;
  d3Shapey[5] = 0;
  d3Shapey[6] = 0;
  d3Shapey[7] = 0;

  d3Shapex2y[0] = -onehalf;
  d3Shapex2y[1] = -onehalf;
  d3Shapex2y[2] = onehalf;
  d3Shapex2y[3] = onehalf;
  d3Shapex2y[4] = 1;
  d3Shapex2y[5] = 0;
  d3Shapex2y[6] = -1;
  d3Shapex2y[7] = 0;

  d3Shapexy2[0] = -onehalf;
  d3Shapexy2[1] = onehalf;
  d3Shapexy2[2] = onehalf;
  d3Shapexy2[3] = -onehalf;
  d3Shapexy2[4] = 0;
  d3Shapexy2[5] = -1;
  d3Shapexy2[6] = 0;
  d3Shapexy2[7] = 1;
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::LocalToGlobalCoord(Scalar *M, Scalar *m, CoordSetT &cs)
{
  // Computes the global X-, Y- and Z-coordinates of a point P.
  //
  // Inputs:  m  = [ξ, η], the local coordinates of P, and
  //          cs = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: M  = [M₀(ξ,η,X,Y,Z), M₁(ξ,η,X,Y,Z), M₂(ξ,η,X,Y,Z)], the global X-, Y- and Z-coordinates of P.

  Scalar Shape[8];
  GetShapeFctVal(Shape,m);

  Scalar X[8], Y[8], Z[8];
  for(int i=0; i<8; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  M[0] = 0.0; M[1] = 0.0; M[2] = 0.0;
  for(int i=0; i<8; i+=4) {
    M[0] += Shape[i]*X[i] + Shape[i+1]*X[i+1] + Shape[i+2]*X[i+2] + Shape[i+3]*X[i+3];
    M[1] += Shape[i]*Y[i] + Shape[i+1]*Y[i+1] + Shape[i+2]*Y[i+2] + Shape[i+3]*Y[i+3];
    M[2] += Shape[i]*Z[i] + Shape[i+1]*Z[i+1] + Shape[i+2]*Z[i+2] + Shape[i+3]*Z[i+3];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dMdx  = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ]
  //          dMdy  = [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η]
  //                  where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[8], dShapey[8];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar X[8], Y[8], Z[8];
  for(int i=0; i<8; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = dMdx[1] = dMdx[2] = 0.0;
  dMdy[0] = dMdy[1] = dMdy[2] = 0.0;
  for(int i=0; i<8; i+=4) {
    dMdx[0] += dShapex[i]*X[i] + dShapex[i+1]*X[i+1] + dShapex[i+2]*X[i+2] + dShapex[i+3]*X[i+3];
    dMdx[1] += dShapex[i]*Y[i] + dShapex[i+1]*Y[i+1] + dShapex[i+2]*Y[i+2] + dShapex[i+3]*Y[i+3];
    dMdx[2] += dShapex[i]*Z[i] + dShapex[i+1]*Z[i+1] + dShapex[i+2]*Z[i+2] + dShapex[i+3]*Z[i+3];

    dMdy[0] += dShapey[i]*X[i] + dShapey[i+1]*X[i+1] + dShapey[i+2]*X[i+2] + dShapey[i+3]*X[i+3];
    dMdy[1] += dShapey[i]*Y[i] + dShapey[i+1]*Y[i+1] + dShapey[i+2]*Y[i+2] + dShapey[i+3]*Y[i+3];
    dMdy[2] += dShapey[i]*Z[i] + dShapey[i+1]*Z[i+1] + dShapey[i+2]*Z[i+2] + dShapey[i+3]*Z[i+3];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs)
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
  Scalar d2Shapex[8], d2Shapey[8], d2Shapexy[8];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar X[8], Y[8], Z[8];
  for(int i=0; i<8; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d2Mdx2[0] = d2Mdx2[1] = d2Mdx2[2] = 0.0;
  d2Mdy2[0] = d2Mdy2[1] = d2Mdy2[2] = 0.0;
  d2Mdxdy[0] = d2Mdxdy[1] = d2Mdxdy[2] = 0.0;
  for(int i=0; i<8; i+=4) {
    d2Mdx2[0] += d2Shapex[i]*X[i] + d2Shapex[i+1]*X[i+1] + d2Shapex[i+2]*X[i+2] + d2Shapex[i+3]*X[i+3];
    d2Mdx2[1] += d2Shapex[i]*Y[i] + d2Shapex[i+1]*Y[i+1] + d2Shapex[i+2]*Y[i+2] + d2Shapex[i+3]*Y[i+3];
    d2Mdx2[2] += d2Shapex[i]*Z[i] + d2Shapex[i+1]*Z[i+1] + d2Shapex[i+2]*Z[i+2] + d2Shapex[i+3]*Z[i+3];

    d2Mdy2[0] += d2Shapey[i]*X[i] + d2Shapey[i+1]*X[i+1] + d2Shapey[i+2]*X[i+2] + d2Shapey[i+3]*X[i+3];
    d2Mdy2[1] += d2Shapey[i]*Y[i] + d2Shapey[i+1]*Y[i+1] + d2Shapey[i+2]*Y[i+2] + d2Shapey[i+3]*Y[i+3];
    d2Mdy2[2] += d2Shapey[i]*Z[i] + d2Shapey[i+1]*Z[i+1] + d2Shapey[i+2]*Z[i+2] + d2Shapey[i+3]*Z[i+3];

    d2Mdxdy[0] += d2Shapexy[i]*X[i] + d2Shapexy[i+1]*X[i+1] + d2Shapexy[i+2]*X[i+2] + d2Shapexy[i+3]*X[i+3];
    d2Mdxdy[1] += d2Shapexy[i]*Y[i] + d2Shapexy[i+1]*Y[i+1] + d2Shapexy[i+2]*Y[i+2] + d2Shapexy[i+3]*Y[i+3];
    d2Mdxdy[2] += d2Shapexy[i]*Z[i] + d2Shapexy[i+1]*Z[i+1] + d2Shapexy[i+2]*Z[i+2] + d2Shapexy[i+3]*Z[i+3];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs)
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
  Scalar d3Shapex[8], d3Shapey[8], d3Shapex2y[8], d3Shapexy2[8];
  Getd3ShapeFct(d3Shapex, d3Shapey, d3Shapex2y, d3Shapexy2, m);

  // Compute ∂³M/∂ξ³, ∂³M/∂η³, ∂³M/∂ξ²∂η and ∂³M/∂ξ∂η²
  Scalar X[8], Y[8], Z[8];
  for(int i=0; i<8; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d3Mdx3[0] = d3Mdx3[1] = d3Mdx3[2] = 0.0;
  d3Mdy3[0] = d3Mdy3[1] = d3Mdy3[2] = 0.0;
  d3Mdx2dy[0] = d3Mdx2dy[1] = d3Mdx2dy[2] = 0.0;
  d3Mdxdy2[0] = d3Mdxdy2[1] = d3Mdxdy2[2] = 0.0;
  for(int i=0; i<8; i+=4) {
/* all zero
    d3Mdx3[0] += d3Shapex[i]*X[i] + d3Shapex[i+1]*X[i+1] + d3Shapex[i+2]*X[i+2] + d3Shapex[i+3]*X[i+3];
    d3Mdx3[1] += d3Shapex[i]*Y[i] + d3Shapex[i+1]*Y[i+1] + d3Shapex[i+2]*Y[i+2] + d3Shapex[i+3]*Y[i+3];
    d3Mdx3[2] += d3Shapex[i]*Z[i] + d3Shapex[i+1]*Z[i+1] + d3Shapex[i+2]*Z[i+2] + d3Shapex[i+3]*Z[i+3];

    d3Mdy3[0] += d3Shapey[i]*X[i] + d3Shapey[i+1]*X[i+1] + d3Shapey[i+2]*X[i+2] + d3Shapey[i+3]*X[i+3];
    d3Mdy3[1] += d3Shapey[i]*Y[i] + d3Shapey[i+1]*Y[i+1] + d3Shapey[i+2]*Y[i+2] + d3Shapey[i+3]*Y[i+3];
    d3Mdy3[2] += d3Shapey[i]*Z[i] + d3Shapey[i+1]*Z[i+1] + d3Shapey[i+2]*Z[i+2] + d3Shapey[i+3]*Z[i+3];
*/
    d3Mdx2dy[0] += d3Shapex2y[i]*X[i] + d3Shapex2y[i+1]*X[i+1] + d3Shapex2y[i+2]*X[i+2] + d3Shapex2y[i+3]*X[i+3];
    d3Mdx2dy[1] += d3Shapex2y[i]*Y[i] + d3Shapex2y[i+1]*Y[i+1] + d3Shapex2y[i+2]*Y[i+2] + d3Shapex2y[i+3]*Y[i+3];
    d3Mdx2dy[2] += d3Shapex2y[i]*Z[i] + d3Shapex2y[i+1]*Z[i+1] + d3Shapex2y[i+2]*Z[i+2] + d3Shapex2y[i+3]*Z[i+3];

    d3Mdxdy2[0] += d3Shapexy2[i]*X[i] + d3Shapexy2[i+1]*X[i+1] + d3Shapexy2[i+2]*X[i+2] + d3Shapexy2[i+3]*X[i+3];
    d3Mdxdy2[1] += d3Shapexy2[i]*Y[i] + d3Shapexy2[i+1]*Y[i+1] + d3Shapexy2[i+2]*Y[i+2] + d3Shapexy2[i+3]*Y[i+3];
    d3Mdxdy2[2] += d3Shapexy2[i]*Z[i] + d3Shapexy2[i+1]*Z[i+1] + d3Shapexy2[i+2]*Z[i+2] + d3Shapexy2[i+3]*Z[i+3];
  }
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceQuad8::GetJacobian(Scalar *m, CoordSetT &cs)
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
FaceQuad8::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  GetShapeFctVal(Shape, m);
  return(GetJacobian(m, cs));
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::GetIsoParamMappingNormalJacobianProduct(Scalar *JNormal, Scalar *m, CoordSetT &cs)
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
FaceQuad8::ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs)
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
FaceQuad8::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs)
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
FaceQuad8::GetIsoParamMappingNormalAndJacobian(Scalar *Normal, Scalar *m, CoordSetT &cs)
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
FaceQuad8::GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.

  // Inputs:  m           = [ξ, η], the local coordinates of P, and
  //          cs          = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormal[0] = [∂n₀/∂X₀, ∂n₀/∂Y₀, ∂n₀/∂Z₀, ∂n₀/∂X₁, ∂n₀/∂Y₁, ∂n₀/∂Z₁, ⋯ , ∂n₀/∂X₇, ∂n₀/∂Y₇, ∂n₀/∂Z₇]
  //          dJNormal[1] = [∂n₁/∂X₀, ∂n₁/∂Y₀, ∂n₁/∂Z₀, ∂n₁/∂X₁, ∂n₁/∂Y₁, ∂n₁/∂Z₁, ⋯ , ∂n₁/∂X₇, ∂n₁/∂Y₇, ∂n₁/∂Z₇]
  //          dJNormal[2] = [∂n₂/∂X₀, ∂n₂/∂Y₀, ∂n₂/∂Z₀, ∂n₂/∂X₁, ∂n₂/∂Y₁, ∂n₂/∂Z₁, ⋯ , ∂n₂/∂X₇, ∂n₂/∂Y₇, ∂n₂/∂Z₇]
  //                        where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                        at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                        where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                        and X = [X₀, X₁, ⋯ , X₇], Y = [Y₀, Y₁, ⋯ , Y₇] and  Z = [Z₀, Z₁, ⋯ , Z₇] are
  //                        the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[8], dShapey[8];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute dJNormal
  for(int i = 0; i < 8; ++i) {
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
FaceQuad8::ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the mixed second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates and the local coordinate
  // of the components of the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m                = [ξ, η], the local coordinates of P, and
  //          cs               = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: ddJNormaldx[][0] = [∂²n₀/∂ξ∂X₀, ∂²n₀/∂ξ∂Y₀, ∂²n₀/∂ξ∂Z₀, ∂²n₀/∂ξ∂X₁, ∂²n₀/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₀/∂ξ∂X₇, ∂²n₀/∂ξ∂Y₇, ∂²n₀/∂ξ∂Z₇]
  //          ddJNormaldx[][1] = [∂²n₁/∂ξ∂X₀, ∂²n₁/∂ξ∂Y₀, ∂²n₁/∂ξ∂Z₀, ∂²n₁/∂ξ∂X₁, ∂²n₁/∂ξ∂Y₁, ∂²n₁/∂ξ∂Z₁, ⋯ , ∂²n₁/∂ξ∂X₇, ∂²n₁/∂ξ∂Y₇, ∂²n₁/∂ξ∂Z₇]
  //          ddJNormaldx[][2] = [∂²n₂/∂ξ∂X₀, ∂²n₂/∂ξ∂Y₀, ∂²n₂/∂ξ∂Z₀, ∂²n₂/∂ξ∂X₁, ∂²n₂/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₂/∂ξ∂X₇, ∂²n₂/∂ξ∂Y₇, ∂²n₂/∂ξ∂Z₇]
  //          ddJNormaldy[][0] = [∂²n₀/∂η∂X₀, ∂²n₀/∂η∂Y₀, ∂²n₀/∂η∂Z₀, ∂²n₀/∂η∂X₁, ∂²n₀/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₀/∂η∂X₇, ∂²n₀/∂η∂Y₇, ∂²n₀/∂η∂Z₇]
  //          ddJNormaldy[][1] = [∂²n₁/∂η∂X₀, ∂²n₁/∂η∂Y₀, ∂²n₁/∂η∂Z₀, ∂²n₁/∂η∂X₁, ∂²n₁/∂η∂Y₁, ∂²n₁/∂η∂Z₁, ⋯ , ∂²n₁/∂η∂X₇, ∂²n₁/∂η∂Y₇, ∂²n₁/∂η∂Z₇]
  //          ddJNormaldy[][2] = [∂²n₂/∂η∂X₀, ∂²n₂/∂η∂Y₀, ∂²n₂/∂η∂Z₀, ∂²n₂/∂η∂X₁, ∂²n₂/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₂/∂η∂X₇, ∂²n₂/∂η∂Y₇, ∂²n₂/∂η∂Z₇]
  //                             where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                             at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                             where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                             and X = [X₀, X₁, ⋯ , X₇], Y = [Y₀, Y₁, ⋯ , Y₇] and  Z = [Z₀, Z₁, ⋯ , Z₇] are
  //                             the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[8], dShapey[8];
  GetdShapeFct(dShapex, dShapey, m);
  
  // Compute shape functions' 2nd derivatives w.r.t. the local coordinates
  Scalar d2Shapex[8], d2Shapey[8], d2Shapexy[8];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ddJNormaldx and ddJNormaldy
  for(int i = 0; i < 8; ++i) {
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

  for(int i = 0; i < 8; ++i) {
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
FaceQuad8::Getd2JNormal(Scalar H[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m    = [ξ, η], the local coordinates of P, and
  //          cs   = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: H[0] = [∂²n₀/∂X₀²,   ∂²n₀/∂X₀∂Y₀, ∂²n₀/∂X₀∂Z₀, ∂²n₀/∂X₀∂X₁, ∂²n₀/∂X₀∂Y₁, ∂²n₀/∂X₀∂Z₁, ⋯ , ∂²n₀/∂X₀∂X₇, ∂²n₀/∂X₀∂Y₇, ∂²n₀/∂X₀∂Z₇,
  //                  ∂²n₀/∂Y₀∂X₀, ∂²n₀/∂Y₀²,   ∂²n₀/∂Y₀∂Z₀, ∂²n₀/∂Y₀∂X₁, ∂²n₀/∂Y₀∂Y₁, ∂²n₀/∂Y₀∂Z₁, ⋯ , ∂²n₀/∂Y₀∂X₇, ∂²n₀/∂Y₀∂Y₇, ∂²n₀/∂Y₀∂Z₇,
  //                                                                         ⋮ 
  //                  ∂²n₀/∂Z₇∂X₀, ∂²n₀/∂Z₇∂Y₀, ∂²n₀/∂Z₇∂Z₀, ∂²n₀/∂Z₇∂X₁, ∂²n₀/∂Z₇∂Y₁, ∂²n₀/∂Z₇∂Z₁, ⋯ , ∂²n₀/∂Z₇∂X₇, ∂²n₀/∂Z₇∂Y₇, ∂²n₀/∂Z₇²]
  //          H[1] = [∂²n₁/∂X₀²,   ∂²n₁/∂X₀∂Y₀, ∂²n₁/∂X₀∂Z₀, ∂²n₁/∂X₀∂X₁, ∂²n₁/∂X₀∂Y₁, ∂²n₁/∂X₀∂Z₁, ⋯ , ∂²n₁/∂X₀∂X₇, ∂²n₁/∂X₀∂Y₇, ∂²n₁/∂X₀∂Z₇,
  //                  ∂²n₁/∂Y₀∂X₀, ∂²n₁/∂Y₀²,   ∂²n₁/∂Y₀∂Z₀, ∂²n₁/∂Y₀∂X₁, ∂²n₁/∂Y₀∂Y₁, ∂²n₁/∂Y₀∂Z₁, ⋯ , ∂²n₁/∂Y₀∂X₇, ∂²n₁/∂Y₀∂Y₇, ∂²n₁/∂Y₀∂Z₇,
  //                                                                         ⋮ 
  //                  ∂²n₁/∂Z₇∂X₀, ∂²n₁/∂Z₇∂Y₀, ∂²n₁/∂Z₇∂Z₀, ∂²n₁/∂Z₇∂X₁, ∂²n₁/∂Z₇∂Y₁, ∂²n₁/∂Z₇∂Z₁, ⋯ , ∂²n₁/∂Z₇∂X₇, ∂²n₁/∂Z₇∂Y₇, ∂²n₁/∂Z₇²]
  //          H[2] = [∂²n₂/∂X₀²,   ∂²n₂/∂X₀∂Y₀, ∂²n₂/∂X₀∂Z₀, ∂²n₂/∂X₀∂X₁, ∂²n₂/∂X₀∂Y₁, ∂²n₂/∂X₀∂Z₁, ⋯ , ∂²n₂/∂X₀∂X₇, ∂²n₂/∂X₀∂Y₇, ∂²n₂/∂X₀∂Z₇,
  //                  ∂²n₂/∂Y₀∂X₀, ∂²n₂/∂Y₀²,   ∂²n₂/∂Y₀∂Z₀, ∂²n₂/∂Y₀∂X₁, ∂²n₂/∂Y₀∂Y₁, ∂²n₂/∂Y₀∂Z₁, ⋯ , ∂²n₂/∂Y₀∂X₇, ∂²n₂/∂Y₀∂Y₇, ∂²n₂/∂Y₀∂Z₇,
  //                                                                         ⋮ 
  //                  ∂²n₂/∂Z₇∂X₀, ∂²n₂/∂Z₇∂Y₀, ∂²n₂/∂Z₇∂Z₀, ∂²n₂/∂Z₇∂X₁, ∂²n₂/∂Z₇∂Y₁, ∂²n₂/∂Z₇∂Z₁, ⋯ , ∂²n₂/∂Z₇∂X₇, ∂²n₂/∂Z₇∂Y₇, ∂²n₂/∂Z₇²]
  //                 where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                 at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                 where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                 and X = [X₀, X₁, ⋯ , X₇], Y = [Y₀, Y₁, ⋯ , Y₇] and  Z = [Z₀, Z₁, ⋯ , Z₇] are
  //                 the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' first partial derivatives w.r.t. the local coordinates
  Scalar dShapex[8], dShapey[8];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute d2JNormal
  for(int i = 0; i < 8; ++ i) {
    for(int j = 0; j < 8; ++j) {
      H[24*(3*j  )+3*i  ][0] = 0;
      H[24*(3*j  )+3*i  ][1] = 0;
      H[24*(3*j  )+3*i  ][2] = 0;

      H[24*(3*j  )+3*i+1][0] = 0;
      H[24*(3*j  )+3*i+1][1] = 0;
      H[24*(3*j  )+3*i+1][2] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];

      H[24*(3*j  )+3*i+2][0] = 0;
      H[24*(3*j  )+3*i+2][1] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[24*(3*j  )+3*i+2][2] = 0;

      H[24*(3*j+1)+3*i  ][0] = 0;
      H[24*(3*j+1)+3*i  ][1] = 0;
      H[24*(3*j+1)+3*i  ][2] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];

      H[24*(3*j+1)+3*i+1][0] = 0;
      H[24*(3*j+1)+3*i+1][1] = 0;
      H[24*(3*j+1)+3*i+1][2] = 0;

      H[24*(3*j+1)+3*i+2][0] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[24*(3*j+1)+3*i+2][1] = 0;
      H[24*(3*j+1)+3*i+2][2] = 0;

      H[24*(3*j+2)+3*i  ][0] = 0;
      H[24*(3*j+2)+3*i  ][1] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[24*(3*j+2)+3*i  ][2] = 0;

      H[24*(3*j+2)+3*i+1][0] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[24*(3*j+2)+3*i+1][1] = 0;
      H[24*(3*j+2)+3*i+1][2] = 0;

      H[24*(3*j+2)+3*i+2][0] = 0;
      H[24*(3*j+2)+3*i+2][1] = 0;
      H[24*(3*j+2)+3*i+2][2] = 0;
    }
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad8::GetUnitNormal(Scalar Normal[3], Scalar *m, CoordSetT &cs)
{
  // Computes the unit vector normal to the surface at a point P.

  GetIsoParamMappingNormalJacobianProduct(Normal, m, cs);

  using std::sqrt;
  Scalar J = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(J != 0.0) {
    Normal[0] /= J; Normal[1] /= J; Normal[2] /= J;
  }
}

#endif
