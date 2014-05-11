#ifndef _THERMBRICK_H_
#define _THERMBRICK_H_

#include <Element.d/Element.h>

class ThermBrick: public Element {

	int nn[8];
public:
	ThermBrick(int*);

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

	Corotator *	getCorotator(CoordSet &cs, double* kel, int, int);

};
#endif

