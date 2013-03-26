#ifndef _HELMQUADGAL_H_
#define _HELMQUADGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGal: public HelmElement, public Element {

	int nn[4];
public:
	HelmQuadGal(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        FullSquareMatrix acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        void getHelmForce(CoordSet&, ComplexVector &, ComplexVector &);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	int getTopNumber();

        PrioInfo examine(int sub, MultiFront *);

        void            addFaces(PolygonSet *pset);

};
#endif

