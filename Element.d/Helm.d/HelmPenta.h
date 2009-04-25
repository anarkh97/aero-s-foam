#ifndef _HELMPENTA_H_
#define _HELMPENTA_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmPenta: public HelmElement, public Element {

	int nn[6];
public:
	HelmPenta(int*);

        Element *clone();

	void renum(int *);

	FullSquareMatrix  stiffness(CoordSet&, double *d, int flg = 1);
	FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	void            addFaces(PolygonSet *pset);

	int		 getTopNumber();

        PrioInfo examine(int sub, MultiFront *mf);
};
#endif

