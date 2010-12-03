#ifndef _QUADPRESSUREBC_H_
#define _QUADPRESSUREBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class QuadPressureBC : public SommerElement
{
    int nnode, nndof, ndime, optele;
    int nn[4];
    double pressure;

  public:
    QuadPressureBC(int *, double);

    int numNodes() { return nnode; }
    int getNode(int nd) { return nn[nd]; }
    int* getNodes() { return nn; }
    int  numDofs();
    int dim() { return ndime; }
    int* dofs(DofSetArray &, int *p=0);
    void markDofs(DofSetArray &);
    void getNormal(CoordSet&, double[3]);

    FullSquareMatrix sommerMatrix(CoordSet&, double *);
    void neumVector(CoordSet&, Vector&, int = 0, GeomState* = 0);
};
#endif

