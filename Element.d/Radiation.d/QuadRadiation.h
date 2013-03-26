#ifndef _QUADRADIATION_H_
#define _QUADRADIATION_H_

#include <Element.d/Element.h>

class QuadRadiation: public virtual Element {

        int nn[4];
public:
        QuadRadiation(int*);

        Element *clone();

        void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg=1);
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

	Corotator* getCorotator(CoordSet &, double*, int, int);

        int              numNodes();
        int*             nodes(int * = 0);
        int              getTopNumber();
        PrioInfo examine(int sub, MultiFront *);
};

#endif

