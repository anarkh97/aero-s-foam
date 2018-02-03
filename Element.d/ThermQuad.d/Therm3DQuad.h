#ifndef _THERM3DQUAD_H_
#define _THERM3DQUAD_H_

#include <Element.d/Element.h>

class Therm3DQuad: public Element {

	int nn[4];
public:
	Therm3DQuad(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg) override;
	FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg) override;
	double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
	int*             dofs(DofSetArray &, int *p=0);
	int              numDofs() const override;

	int             numNodes() const override;
	int*             nodes(int * = 0) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;
	void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
	                       int hflInd);
	Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};

#endif
