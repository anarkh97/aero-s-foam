#ifndef _HELMTRI3GLS_H_
#define _HELMTRI3GLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmTri3Gls: public HelmElement, public Element {

	int nn[3];
	mutable double coef; // TODO Get rid of this variable.
public:
	HelmTri3Gls(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg=1) const;
	FullSquareMatrix  acousticm(CoordSet& cs, double *d);
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const;

	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p=0) override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void            addFaces(PolygonSet *pset);
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;

	virtual double getHelmCoef() { return coef; }

};
#endif

