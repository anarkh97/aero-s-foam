#ifndef _THERM3NOSHELL_H_
#define _THERM3NOSHELL_H_

#include	<Element.d/Element.h>

class GeomState;

class Therm3NoShell : public Element {

	int nn[3];
public:
	explicit Therm3NoShell(int*);

	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;

	void computeTemp(CoordSet&cs, State &state, double gp[2], double*res) override;

	void getFlFlux(double gp[2], double *flF, double *resF) override;
	void getGravityForce(CoordSet& cs, double *, Vector &force, int, GeomState *) override;
	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};
#endif

