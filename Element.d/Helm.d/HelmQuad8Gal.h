#ifndef _HELMQUAD8GAL_H_
#define _HELMQUAD8GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuad8Gal: public HelmElement, public Element {

	int nn[8];
public:
	HelmQuad8Gal(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix acousticm(CoordSet&, double *d);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double           getMass(CoordSet&);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override { return(132); }
	void		addFaces(PolygonSet *pset);

};
#endif

