#ifndef _MATRIXELEMENT_H_
#define _MATRIXELEMENT_H_

#include <Element.d/Element.h>

class DofSet;
template <class Scalar> class GenAssembledFullM;

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
	~MatrixElement() override;

    Element* clone() override;

    void setDofs(DofSet *d);
    void setStiffness(GenAssembledFullM<double> *k);
    void setStiffness(GenAssembledFullM<complex<double> > *k);

    void renum(int *table) override;
    void renum(EleRenumMap& m) override;
    void markDofs(DofSetArray &) override;
    int* dofs(DofSetArray &, int *p) override;
    int* nodes(int *) const override;

    FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
    FullSquareMatrix imagStiffness(CoordSet&, double *kel, int flg=1);

    int  numDofs() const { return ndofs; }
    int  numNodes() const { return nnodes; }

};
#endif

