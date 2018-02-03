#ifndef _BRICKCONTACT_H_
#define _BRICKCONTACT_H_

#include <Element.d/Element.h>

class BrickContact: public Element {

	int nn[8];
public:
	BrickContact(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const;
        FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const;
	double           getMass(const CoordSet& cs) const;


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;

};
#endif

