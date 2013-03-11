#ifndef _QUADPRESSUREBC_H_
#define _QUADPRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Quad4LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Quad4LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Quad4LagrangePolynomialShapeFunction, GaussLegendre2d>;

class QuadPressureBC : public PressureElement<Quad4LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    QuadPressureBC(int* _nn, double _pressure, bool); 

  protected:
    double pressure;
    void getConstants(CoordSet& cs, Eigen::Array<double,13,1>&, Eigen::Array<int,1,1>&);
};

#else
#include <Element.d/Sommerfeld.d/SommerElement.h>

class QuadPressureBC : public SommerElement
{
    int nnode, nndof, ndime, optele;
    int nn[4];
    double pressure;
    bool ConwepOnOff;

  public:
    QuadPressureBC(int *, double, bool _ConwepOnOff);

    int numNodes() { return nnode; }
    int getNode(int nd) { return nn[nd]; }
    int* getNodes() { return nn; }
    int  numDofs();
    int dim() { return ndime; }
    int* dofs(DofSetArray &, int *p=0);
    void markDofs(DofSetArray &);
    void getNormal(CoordSet&, double[3]);

    FullSquareMatrix sommerMatrix(CoordSet&, double *);
    void neumVector(CoordSet&, Vector&, int = 0, GeomState* = 0, double time = 0.0);

    int findAndSetEle(CoordSet& cs,Elemset &eset,
        Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
        int it = 0) { return 0; } // normals will never be flipped. TODO reconsider this
};

#endif
#endif
