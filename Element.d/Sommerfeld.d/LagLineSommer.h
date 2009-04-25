#ifndef _LAGLINESOMMER_H_
#define _LAGLINESOMMER_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class LagLineSommer: public SommerElement {

	int *nn;
        int order;
public:
	LagLineSommer(int, int*, Element *_el = 0);

        int numNodes() {return order;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return order; }
        int dim() { return 2; }
        int* dofs(DofSetArray &, int *p=0);
        virtual LagLineSommer* clone();

        void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double);

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return exp(ComplexD(0.0,M_PI/4.0))/sqrt(8.0*M_PI*k); }

	void getNormal(CoordSet&, double [3]);
};
#endif

