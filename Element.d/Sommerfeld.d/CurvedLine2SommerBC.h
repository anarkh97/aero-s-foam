#ifndef _CURVEDLINE2SOMMERBC_H_
#define _CURVEDLINE2SOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class CurvedLine2SommerBC: public SommerElement {

	int nn[3];
public:
	CurvedLine2SommerBC(int, int, int, Element *_el = 0);

        int numNodes() {return 3;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return 3; }
        int dim() { return 2; }
        int* dofs(DofSetArray &, int *p=0);
        virtual CurvedLine2SommerBC* clone();

        void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int pflag=0);

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return exp(ComplexD(0.0,M_PI/4.0))/sqrt(8.0*M_PI*k); }

	void getNormal(CoordSet&, double [3]);
};
#endif

