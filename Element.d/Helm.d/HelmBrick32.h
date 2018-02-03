// ---------------------------------------------------------------------
// HB - 05-24-05
// ---------------------------------------------------------------------
// 32 nodes brick element
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#ifndef _HELMBRICK32_H_
#define _HELMBRICK32_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmBrick32: public HelmElement, public Element {

	int nn[32];
public:
	HelmBrick32(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(CoordSet&, double *d, int flg = 1);
	FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;

	void            addFaces(PolygonSet *pset);

        int getTopNumber() override;
        int numTopNodes();

        PrioInfo examine(int sub, MultiFront *mf) override;
};
#endif

