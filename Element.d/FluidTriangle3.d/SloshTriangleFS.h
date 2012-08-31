#ifndef _SLOSHTRIANGLEFS_H_
#define _SLOSHTRIANGLEFS_H_

#include <Element.d/Element.h>

class SloshTriangleFS: public Element {

	int nn[3];
public:
	SloshTriangleFS(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix  stiffness(CoordSet& cs, double *d, int flg = 1);
        FullSquareMatrix  massMatrix(CoordSet& cs, double *mel, int cmflg=1);
        double            getMass(CoordSet&);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int              getTopNumber();

	PrioInfo examine(int sub, MultiFront *)   {
           fprintf(stderr,"SloshTriangleFS.h: PrioInfo examine is commented in Dec.d/ElemFSCheck.C"); return *(new PrioInfo);
        };
        //void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
        //void getFlFlux(double gp[2], double *flF, double *resF);
        //void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }

};
#endif

