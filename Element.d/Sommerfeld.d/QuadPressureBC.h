#ifndef _QUADPRESSUREBC_H_
#define _QUADPRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Function.d/ExternalForce.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Quad4LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Quad4LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Quad4LagrangePolynomialShapeFunction, GaussLegendre2d>;

class QuadPressureBC : public PressureElement<Quad4LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    QuadPressureBC(int* _nn, PressureBCond* _pbc); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,21,1>&, Eigen::Array<int,2,1>&);
};

#else
#include <Element.d/Sommerfeld.d/SommerElement.h>

class QuadPressureBC : public SommerElement
{
    int nnode, nndof, ndime, optele;
    int nn[4];
    PressureBCond *pbc;

  public:
    QuadPressureBC(int *, PressureBCond *);

    int numNodes() { return nnode; }
    int getNode(int nd) { return nn[nd]; }
    int* getNodes() { return nn; }
    int  numDofs();
    int dim() { return ndime; }
    int* dofs(DofSetArray &, int *p=0);
    void markDofs(DofSetArray &);
    void getNormal(CoordSet&, double[3]);

    FullSquareMatrix sommerMatrix(CoordSet&, double *);
    PressureBCond* getPressure() { return pbc; }
    void neumVector(CoordSet&, Vector&, int = 0, GeomState* = 0, double t = 0);

    int findAndSetEle(CoordSet& cs,Elemset &eset,
        Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
        int it = 0) { return 0; } // normals will never be flipped. TODO reconsider this
};

#endif
#endif
