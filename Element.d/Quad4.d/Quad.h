#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>

class Quad: public Element {

        int numnod; // number of nodes in the element
	int *nn;
public:
	Quad(int , int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, 
	                                 int gravflg, GeomState *gs);
        virtual void     getVonMises (Vector &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp,
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	int getTopNumber();

};
#endif

