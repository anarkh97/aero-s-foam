// ---------------------------------------------------------------------
// HB - 05-24-05
// ---------------------------------------------------------------------
// 26 nodes wedge element
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#ifndef _HELMPENTA26_H_
#define _HELMPENTA26_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmPenta26: public HelmElement, public Element {

	int nn[26];
public:
	HelmPenta26(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg = 1) const;
	FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
	double           getMass(const CoordSet& cs) const;

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

