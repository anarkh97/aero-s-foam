#ifndef _HELMISOPARAMTETRA_H_
#define _HELMISOPARAMTETRA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamTetra: public HelmElement, public Element {

	int order;
	int *nn;
	HelmIsoParamTetra(const HelmIsoParamTetra& e);

public:
	HelmIsoParamTetra(int,int*);

	FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
	FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
	FullSquareMatrixC stiffness(CoordSet&, complex<double> *d);
	FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d);
	double getMass(CoordSet&);

	Element *clone() override;
	void renum(int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
	int getTopNumber() override {return 196;}
	int numTopNodes() {return (order*(order+1)*(order+2))/6;}
	int* dofs(DofSetArray &, int *p) override;
	int numDofs() const { return (order*(order+1)*(order+2))/6; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;
	void addFaces(PolygonSet *pset) override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	int nDecFaces() { return 4;}
	int getDecFace(int iFace, int *fn);

};
#endif
