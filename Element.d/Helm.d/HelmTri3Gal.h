#ifndef _HELMTRI3GAL_H_
#define _HELMTRI3GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix22[2][2];

class HelmTri3Gal: public HelmElement, public Element {

	int nn[3];
public:
	HelmTri3Gal(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg=1) const;
        FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;

        void getHelmForce(CoordSet& cs, ComplexVector &vc, ComplexVector &force);

        void computedxdxi(CoordSet &cs, int nint, double (*derivatives)[3][2],
                          Matrix22 *dxdxi, double *det);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
                            double kappa, double *waveDir);

        double           getMass(const CoordSet&) const;

	void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p=0) const override;
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;

	void            addFaces(PolygonSet *pset);
		PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;

};

#endif

