#ifndef _HELMLAGQUADGAL_H_
#define _HELMLAGQUADGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmLagQuadGal: public HelmElement, public Element {

	int order;
	int *nn;
	void shapeFunctions(double xi, double eta, double *N);
	HelmLagQuadGal(const HelmLagQuadGal& e);

public:
	HelmLagQuadGal(int,int*);

	FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
	FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
	FullSquareMatrix acousticm(CoordSet&, double *d);
	void wErrors(CoordSet&,
	             double *l2e, double *h1e, double *l2, double *h1,
	             ComplexD *u, double kappa, double *waveDir);
	double           getMass(CoordSet&);
	void edgeShapeFunctions(int n1, int n2, int *ng,
	                        double **gw, double **N);

	Element *clone() override;
	void renum(int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	int numDofs() const override { return order*order; }
	int numNodes() const override { return order*order; }
	int * nodes(int *) const override;
	void addFaces(PolygonSet *pset) override;
	int getTopNumber() override {return 163;}

	PrioInfo examine(int sub, MultiFront *mf) override;
};
#endif

