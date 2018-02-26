#ifndef _THERMSOPARAMHEXA_H_
#define _THERMSOPARAMHEXA_H_

#include <Element.d/Element.h>

class ThermIsoParamHexa: public Element {

	int order;
	int *nn;
	ThermIsoParamHexa(const ThermIsoParamHexa& e);

public:
	ThermIsoParamHexa(int,int*);

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;

	Element *clone() override;
	void renum(int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int getTopNumber() override {return 195;}
	int numTopNodes() const override {return order*order*order;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return order*order*order; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;

	PrioInfo examine(int sub, MultiFront *mf) override;

	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};
#endif
