#ifndef _HELMBRICKGLS_H_
#define _HELMBRICKGLS_H_

#include <cmath>
#include <Element.d/Helm.d/HelmElement.h>

class HelmBrickGLS: public HelmElement, public Element {

	int nn[8];
        double coef;
public:
	HelmBrickGLS(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(CoordSet&, double *d, int flg=1);
	FullSquareMatrix  acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

	void            addFaces(PolygonSet *pset);

	int getTopNumber() override;

        virtual double helmCoef() { return coef; }

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() { return 6;}
        int getDecFace(int iFace, int *fn);
};
#endif

