#ifndef _HELMISOPARAMQUAD_H_
#define _HELMISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamQuad: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamQuad(const HelmIsoParamQuad& e);

public:
	HelmIsoParamQuad(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
        FullSquareMatrixC stiffness(const CoordSet&, complex<double> *d) const;
        FullSquareMatrixC massMatrix(const CoordSet&, complex<double> *d) const;
        double getMass(const CoordSet&) const;

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
//	int getTopNumber() override {return 195;}
	int numTopNodes() {return order*order;}
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const override { return order*order; }
        int numNodes() const;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;

};
#endif
