#ifndef _TRIANGLECONVEC_H_
#define _TRIANGLECONVEC_H_

#include <Element.d/Element.h>

class TriangleConvec: public Element {

	int nn[4];
public:
	TriangleConvec(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg=1);
        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int		 getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
};

#endif
