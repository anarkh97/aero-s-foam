#ifndef _HELMQUADGAL_H_
#define _HELMQUADGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGal: public HelmElement, public Element {

	int nn[4];
public:
	HelmQuadGal(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg = 1) const;
        FullSquareMatrix acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const;
        void getHelmForce(CoordSet&, ComplexVector &, ComplexVector &);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

	int getTopNumber() override;

        PrioInfo examine(int sub, MultiFront *) override;

        void            addFaces(PolygonSet *pset);

};
#endif

