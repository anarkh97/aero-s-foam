#ifndef _HELMPENTA_H_
#define _HELMPENTA_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmPenta: public HelmElement, public Element {

	int nn[6];
public:
	HelmPenta(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg = 1) const;
	FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
	double           getMass(const CoordSet& cs) const;

	void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p=0) const override;
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

	void            addFaces(PolygonSet *pset);

	int getTopNumber() override;

        PrioInfo examine(int sub, MultiFront *mf) override;
};
#endif

