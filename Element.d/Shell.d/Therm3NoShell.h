#ifndef _THERM3NOSHELL_H_
#define _THERM3NOSHELL_H_

#include	<Element.d/Element.h>

class GeomState;

class Therm3NoShell : public Element {

	int nn[3];
public:
	Therm3NoShell(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int              getTopNumber();
        bool hasRot() {return true;} // DEC
	PrioInfo examine(int sub, MultiFront *);

        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
                        
	void getFlFlux(double gp[2], double *flF, double *resF);

        //bool isShell() { return true; }

        void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }

};
#endif

