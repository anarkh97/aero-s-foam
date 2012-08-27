#ifndef _HELMQUADGLS_H_
#define _HELMQUADGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGls: public HelmElement, public Element {

	int nn[4];
        double coef;
public:
	HelmQuadGls(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix  stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix  acousticm(CoordSet& cs, double *d);
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg=1);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();

	PrioInfo examine(int sub, MultiFront *);

	void            addFaces(PolygonSet *pset);

        virtual double helmCoef() { return coef; }
};
#endif

