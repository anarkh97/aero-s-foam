#ifndef _HELMAXITRI_H_
#define _HELMAXITRI_H_

#include <HelmAxi.d/AxiHElem.h>

class PolygonSet;

class HelmAxiTri: public AxiHElement {

	int nn[3];
public:
	HelmAxiTri(int*);
	//AxiHElement *clone();
	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
        FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        FullSquareMatrix stiffteta(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&,double *d,int cmflg = 1);
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

