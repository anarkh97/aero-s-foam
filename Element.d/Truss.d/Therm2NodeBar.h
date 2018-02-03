#ifndef _THERM2NODEBAR_H_
#define _THERM2NODEBAR_H_

#include <Element.d/Element.h>

class Therm2NodeBar: public Element {

	int nn[2];
public:

	explicit Therm2NodeBar(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet&,double *kel, int flg) override;
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg) override;
	double           getMass(CoordSet&) override;

	void             markDofs(DofSetArray &) override;
	int*             dofs(DofSetArray &, int *p) override;
	int              numDofs() const override;

	int             numNodes() const override;
	int*             nodes(int * ) const override;

	int 		getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;
	void computeTemp(CoordSet&cs, State &state, double gp[2], double*res) override;
	void getFlFlux(double gp[2], double *flF, double *resF) override;
	void trussHeatFluxes(double &elheatflux, CoordSet&, Vector &elTemp,
	                     int hflIndex) override;
	void getGravityForce(CoordSet& cs, double *, Vector &force, int, GeomState *) override;

	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};
#endif
