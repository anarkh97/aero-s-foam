#ifndef _THERM2NODEBAR_H_
#define _THERM2NODEBAR_H_

#include <Element.d/Element.h>

class Therm2NodeBar: public Element {

        int nn[2];
public:

	Therm2NodeBar(int*);

        Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	int 		getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
        void getFlFlux(double gp[2], double *flF, double *resF);
        void trussHeatFluxes(double &elheatflux, CoordSet&, Vector &elTemp, 
                                         int hflIndex);
};
#endif
