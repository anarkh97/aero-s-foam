#ifndef _MATRIXELEMENT_H_
#define _MATRIXELEMENT_H_

#include <Element.d/Element.h>

class DofSet;

class MatrixElement : public Element
{
    int nnodes;  // number of nodes
    int *nn;  // node numbers
    int ndofs; // number of dofs
    DofSet *alldofs;
    GenAssembledFullM<double> *k_real; // stiffness matrix
    GenAssembledFullM<complex<double> > *k_complex;
    int *renumTable;

  public:
    MatrixElement(int _nnodes, int *_nn);
    virtual ~MatrixElement();

    Element* clone();

    void setDofs(DofSet *d);
    void setStiffness(GenAssembledFullM<double> *k);
    void setStiffness(GenAssembledFullM<complex<double> > *k);

    void renum(int *table);
    void renum(EleRenumMap& m);
    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int* nodes(int * = 0);

    FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix imagStiffness(CoordSet&, double *kel, int flg=1);

    int  numDofs() { return ndofs; }
    int  numNodes() { return nnodes; }

};
#endif

