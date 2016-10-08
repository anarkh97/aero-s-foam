#ifndef _LINSPRING_H_
#define _LINSPRING_H_

#include <Element.d/Element.h>

class LinSpring : public Element {

        int nn[1];
public:

	LinSpring(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	bool isSafe()    {return true;}
	bool isSpring()  {return true;}
	int getTopNumber(){return 111;}
	PrioInfo examine(int sub, MultiFront *);
        Corotator *getCorotator(CoordSet &, double*, int, int);
};
#endif
