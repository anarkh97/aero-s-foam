#ifndef _TRIANGLERADIATION_H_
#define _TRIANGLERADIATION_H_

#include <Element.d/Element.h>

class TriangleRadiation: public virtual Element {

        int nn[3];
public:
        TriangleRadiation(int*);
        ~TriangleRadiation();

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

        bool isRadiationElement() { return true; }

        void computeTemp(CoordSet&, State &, double[2], double*); 
        void getFlFlux(double[2], double *, double *);
};

#endif
