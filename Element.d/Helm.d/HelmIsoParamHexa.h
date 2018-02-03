#ifndef _HELMISOPARAMHEXA_H_
#define _HELMISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamHexa: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamHexa(const HelmIsoParamHexa& e);

public:
	HelmIsoParamHexa(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const;
        FullSquareMatrixC stiffness(const CoordSet&, complex<double> *d) const;
        FullSquareMatrixC massMatrix(const CoordSet&, complex<double> *d) const;
        double getMass(const CoordSet&) const;

	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
	int getTopNumber() override {return 195;}
	int numTopNodes() {return order*order*order;}
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const { return order*order*order; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;
	void addFaces(PolygonSet *pset) override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() { return 6;}
        int getDecFace(int iFace, int *fn);
};
#endif
