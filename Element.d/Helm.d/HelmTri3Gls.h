#ifndef _HELMTRI3GLS_H_
#define _HELMTRI3GLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmTri3Gls: public HelmElement, public Element {

	int nn[3];
        double coef;
public:
	HelmTri3Gls(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix  stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix  acousticm(CoordSet& cs, double *d);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	void            addFaces(PolygonSet *pset);
	int getTopNumber();
	PrioInfo examine(int sub, MultiFront *);

        virtual double getHelmCoef() { return coef; }

};
#endif

