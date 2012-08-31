#ifndef _TRIANGLEPRESSUREBC_H_
#define _TRIANGLEPRESSUREBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class TrianglePressureBC : public SommerElement
{
    int nnode, nndof, ndime, optele;
    int nn[3];
    double pressure;

  public:
    TrianglePressureBC(int *, double);

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

    int findAndSetEle(CoordSet& cs,Elemset &eset,
        Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
        int it = 0) { return 0; } // normals will never be flipped. TODO reconsider this
};
#endif

