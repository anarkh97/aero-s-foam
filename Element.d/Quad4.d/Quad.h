#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>

class Quad: public Element {

        int numnod; // number of nodes in the element
	int *nn;
public:
	Quad(int , int*);

	Element *clone() override;

	void renum(const int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg = 1) const override;
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
        double           getMass(const CoordSet&) const override;

        virtual void     getVonMises (Vector &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp,
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);


	void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p) const override;
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

	int getTopNumber() override;

        int getMassType() const override { return 0; } // lumped only
};
#endif

