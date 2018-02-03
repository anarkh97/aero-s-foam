#ifndef _HELMQUADGLS_H_
#define _HELMQUADGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGls: public HelmElement, public Element {

	int nn[4];
        double coef;
public:
	HelmQuadGls(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix  stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix  acousticm(CoordSet& cs, double *d);
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg) override;

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;

	PrioInfo examine(int sub, MultiFront *) override;

	void            addFaces(PolygonSet *pset);

        virtual double helmCoef() { return coef; }
};
#endif

