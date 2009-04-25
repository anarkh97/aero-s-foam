#ifndef _ROTNSPRLINK_H_
#define _ROTNSPRLINK_H_

#include <Element.d/Element.h>

class RotnSprlink : public Element {

        int nn[2];
public:

	RotnSprlink(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet& cs, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

        Corotator *getCorotator(CoordSet &, double*, int, int);
	int getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
        bool isSpring() { return true; }
	
};
#endif
