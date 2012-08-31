#ifndef _HELMTRI6GAL_H_
#define _HELMTRI6GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix22[2][2];

class HelmTri6Gal: public HelmElement, public Element {

	int nn[6];
public:
	HelmTri6Gal(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	void            addFaces(PolygonSet *pset);

        void computedxdxi(CoordSet &cs, int nint, double (*derivatives)[6][2],
                          Matrix22 *dxdxi, double *det);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
                            double kappa, double *waveDir);
	int getTopNumber() {return 138;}
};
#endif

