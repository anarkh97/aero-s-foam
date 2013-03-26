#ifndef _THERMQUADGAL_H_
#define _THERMQUADGAL_H_

#include <Element.d/Element.h>

class ThermQuadGal: public Element {

	int nn[4];
public:
	ThermQuadGal(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flag=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        int		 getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
        void getFlFlux(double gp[2], double *flF, double *resF);
        void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp, 
                               int hflInd);
        void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }
};
#endif

