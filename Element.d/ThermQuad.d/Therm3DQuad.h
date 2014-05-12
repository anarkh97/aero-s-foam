#ifndef _THERM3DQUAD_H_
#define _THERM3DQUAD_H_

#include <Element.d/Element.h>

class Therm3DQuad: public Element {

	int nn[4];
public:
	Therm3DQuad(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

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
        void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
                               int hflInd);
        Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};

#endif
