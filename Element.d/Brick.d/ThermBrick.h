#ifndef _THERMBRICK_H_
#define _THERMBRICK_H_

#include <Element.d/Element.h>

class ThermBrick: public Element {

	int nn[8];
public:
	ThermBrick(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const;
        FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const;
	double           getMass(const CoordSet& cs) const;


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;

	Corotator *	getCorotator(CoordSet &cs, double* kel, int, int);

};
#endif

