#ifndef _HELMTRI3GAL_H_
#define _HELMTRI3GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix22[2][2];

class HelmTri3Gal: public HelmElement, public Element {

	int nn[3];
public:
	HelmTri3Gal(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix  stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

        void getHelmForce(CoordSet& cs, ComplexVector &vc, ComplexVector &force);

        void computedxdxi(CoordSet &cs, int nint, double (*derivatives)[3][2],
                          Matrix22 *dxdxi, double *det);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
                            double kappa, double *waveDir);

        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	void            addFaces(PolygonSet *pset);
		PrioInfo examine(int sub, MultiFront *);
	int getTopNumber();

};

#endif

