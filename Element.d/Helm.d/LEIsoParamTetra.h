#ifndef _LEISOPARAMTETRA_H_
#define _LEISOPARAMTETRA_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamTetra: public Element {

        int order;
	int *nn;
        LEIsoParamTetra(const LEIsoParamTetra& e);

public:
	LEIsoParamTetra(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
        void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega) override;
        double  getMass(const CoordSet& cs) const override;

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p) const override;
        int numDofs() const override { return (order*(order+1)*(order+2))/2; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() const override { return 4;}
        int getDecFace(int iFace, int *fn) override;

};
#endif
