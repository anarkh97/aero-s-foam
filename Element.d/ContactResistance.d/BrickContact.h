#ifndef _BRICKCONTACT_H_
#define _BRICKCONTACT_H_

#include <Element.d/Element.h>

class BrickContact: public Element {

	int nn[8];
public:
	BrickContact(int*);

        Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int		 getTopNumber();
	PrioInfo examine(int sub, MultiFront *);

};
#endif

