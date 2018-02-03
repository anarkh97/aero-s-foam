#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>

class Quad: public Element {

        int numnod; // number of nodes in the element
	int *nn;
public:
	Quad(int , int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        virtual void     getVonMises (Vector &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp,
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

	int getTopNumber() override;

        int getMassType() { return 0; } // lumped only
};
#endif

