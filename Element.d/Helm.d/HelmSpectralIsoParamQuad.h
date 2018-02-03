#ifndef _HELMSPECTRALISOPARAMQUAD_H_
#define _HELMSPECTRALISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmSpectralIsoParamQuad: public HelmElement, public Element {

        int order;
	int *nn;
        HelmSpectralIsoParamQuad(const HelmSpectralIsoParamQuad& e);

public:
	HelmSpectralIsoParamQuad(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const;
        double getMass(const CoordSet&) const;

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
	int getTopNumber() override {return 163;}
	int numTopNodes() {return order*order;}
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const { return order*order; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;

};
#endif
