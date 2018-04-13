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

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	FullSquareMatrixC stiffness(const CoordSet&, complex<double> *d) const override;
	FullSquareMatrixC massMatrix(const CoordSet&, complex<double> *d) const override;
	double  getMass(const CoordSet& cs) const override;

	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int getTopNumber() override {return 196;}
	int numTopNodes() const override {return (order*(order+1)*(order+2))/6;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return (order*(order+1)*(order+2))/6; }
	int numNodes() const override;
	int* nodes(int *) const override;
	void addFaces(PolygonSet *pset) override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	int nDecFaces() const override { return 4;}
	int getDecFace(int iFace, int *fn) override;

};
#endif
