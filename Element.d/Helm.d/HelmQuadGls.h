#ifndef _HELMQUADGLS_H_
#define _HELMQUADGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGls: public HelmElement, public Element {

	int nn[4];
	mutable double coef; // TODO Get rid of this variable.
public:
	HelmQuadGls(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg=1) const override;
	FullSquareMatrix  acousticm(CoordSet& cs, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p=0) override;
	int numDofs() const override;

	int numNodes() const override;
	int *nodes(int *) const override;
	int getTopNumber() override;

	PrioInfo examine(int sub, MultiFront *) override;

	void            addFaces(PolygonSet *pset);

	virtual double helmCoef() { return coef; }
};
#endif

