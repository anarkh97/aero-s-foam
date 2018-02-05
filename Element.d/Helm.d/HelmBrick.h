#ifndef _HELMBRICK_H_
#define _HELMBRICK_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmBrick: public HelmElement, public Element {

	int nn[8];
public:
	HelmBrick(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg = 1) const;
	FullSquareMatrix acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
	double getMass(const CoordSet& cs) const;

	void markDofs(DofSetArray &) override;
        int* dofs(DofSetArray &, int *p) override;
         int numDofs() const override;

        int numNodes() const override;
        int* nodes(int * = 0) const override;

	void addFaces(PolygonSet *pset) override;

	int getTopNumber() override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() const override { return 6;}
        int getDecFace(int iFace, int *fn);
};
#endif

