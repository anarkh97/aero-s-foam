#ifndef _CTCVIRTUALELT_H_
#define _CTCVIRTUALELT_H_

#include <Element.d/Element.h>

class CtcVirtualElt : public Element {

        int nn[2];

public:

	CtcVirtualElt(int,int*);

        Element *clone();

	void renum(int *);
	void markDofs(DofSetArray &);
	int* dofs(DofSetArray &, int *p=0);
	int  numDofs();
	
	int  numNodes();
	int* nodes(int * = 0);
	int  getTopNumber();


	FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
};

#endif
