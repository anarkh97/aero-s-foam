#ifndef _THERMQUADGAL_H_
#define _THERMQUADGAL_H_

#include <Element.d/Element.h>

class ThermQuadGal: public Element {

	int nn[4];
public:
	explicit ThermQuadGal(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet&, double *d, int flag) override;
	FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg) override;
	double           getMass(CoordSet&) override;

	void             markDofs(DofSetArray &) override;
	int*             dofs(DofSetArray &, int *p) override;
	int              numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	int	getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;
	void computeTemp(CoordSet&cs, State &state, double gp[2], double*res) override;
	void getFlFlux(double gp[2], double *flF, double *resF) override;
	void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
	                       int hflInd) override;
	Corotator * getCorotator(CoordSet &, double*, int, int) override { return nullptr; }
};
#endif

