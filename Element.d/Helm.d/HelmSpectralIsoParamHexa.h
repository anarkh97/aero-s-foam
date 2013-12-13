#ifndef _HELMSPECTRALISOPARAMHEXA_H_
#define _HELMSPECTRALISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmSpectralIsoParamHexa: public HelmElement, public Element {

        int order;
	int *nn;
        HelmSpectralIsoParamHexa(const HelmSpectralIsoParamHexa& e);

public:
	HelmSpectralIsoParamHexa(int,int*);
        ~HelmSpectralIsoParamHexa() { delete [] nn; }

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
	int getTopNumber() {return 195;}
	int numTopNodes() {return order*order*order;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return order*order*order; }
        int numNodes();
        int* nodes(int * = 0);
	void addFaces(PolygonSet *pset);

        PrioInfo examine(int sub, MultiFront *mf);

};
#endif
