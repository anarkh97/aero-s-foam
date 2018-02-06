#ifndef _THERMISOPARAMTETRA_H_
#define _THERMISOPARAMTETRA_H_

#include <Element.d/Element.h>

class ThermIsoParamTetra: public Element {

        int order;
	int *nn;
        ThermIsoParamTetra(const ThermIsoParamTetra& e);

public:
	ThermIsoParamTetra(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
        double getMass(const CoordSet&) const;

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const override { return (order*(order+1)*(order+2))/6; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int getTopNumber() override;

        Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};
#endif
