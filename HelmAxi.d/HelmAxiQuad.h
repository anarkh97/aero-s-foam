#ifndef _HELMAXIQUAD_H_
#define _HELMAXIQUAD_H_

#include <HelmAxi.d/AxiHElem.h>

class PolygonSet;

class HelmAxiQuad: public AxiHElement {

	int nn[4];
public:
	HelmAxiQuad(int*);
	//AxiHElement *clone();
	Element *clone();
	void renum(int *);
        FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        FullSquareMatrix stiffteta(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg = 1);
	void    markDofs(DofSetArray &);
        int*    dofs(DofSetArray &, int *p=0);
        int     numDofs();
        int     numNodes();
        int*    nodes(int * = 0);
	void	addFaces(PolygonSet *pset);
        void    buildMesh3D(int &elemNum, FILE *outF,
                            int nodeInc, int numSlices);
        void    buildMesh2D(int &elemNum, FILE *outF,
                            int nodeInc, int numSlices);
	PrioInfo examine(int sub, MultiFront *);
};
#endif

