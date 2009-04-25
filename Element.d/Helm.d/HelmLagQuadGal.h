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

	Element *clone();
	void renum(int *);
	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs() { return order*order; }
        int              numNodes() { return order*order; }
        int*             nodes(int * = 0);
	void		addFaces(PolygonSet *pset);
	int getTopNumber() {return 163;}
	
	double weight() {
	  return order;
	}
	
	double trueWeight() {
	  return weight();
	}

	PrioInfo examine(int sub, MultiFront *mf);
};
#endif

