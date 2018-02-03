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

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega);
        double getMass(CoordSet&);

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const { return (order*(order+1)*(order+2))/2; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() { return 4;}
        int getDecFace(int iFace, int *fn);

};
#endif
