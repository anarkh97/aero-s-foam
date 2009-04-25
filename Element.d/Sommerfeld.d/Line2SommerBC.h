#ifndef _LINE2SOMMERBC_H_
#define _LINE2SOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class Line2SommerBC: public SommerElement {

	int nn[3];
public:
	Line2SommerBC(int, int, int, Element *_el = 0, int etype = 2);

        int numNodes() {return 3;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return 3; }
        int dim() { return 2; }
        int* dofs(DofSetArray &, int *p=0);
        virtual Line2SommerBC* clone();

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return exp(ComplexD(0.0,M_PI/4.0))/sqrt(8.0*M_PI*k); }

	void getNormal(CoordSet&, double [3]);
};
#endif

