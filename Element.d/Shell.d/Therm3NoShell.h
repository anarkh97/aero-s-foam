#ifndef _THERM3NOSHELL_H_
#define _THERM3NOSHELL_H_

#include	<Element.d/Element.h>

class GeomState;

class Therm3NoShell : public Element {

	int nn[3];
public:
	Therm3NoShell(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;

        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
                        
	void getFlFlux(double gp[2], double *flF, double *resF);
        void getGravityForce(CoordSet& cs, double *, Vector &force, int, GeomState *);

        Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};
#endif

