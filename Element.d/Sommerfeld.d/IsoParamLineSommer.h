#ifndef _ISOPARAMLINESOMMER_H_
#define _ISOPARAMLINESOMMER_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamLineSommer: public SommerElement {

	int *nn;
        int order;
public:
	IsoParamLineSommer(int, int*, Element *_el = 0);

        int numNodes() {return order;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return order; }
        int numWetDofs() { return 3*order; }
        int dim() { return 2; }
        int* dofs(DofSetArray &, int *p=0);
        virtual IsoParamLineSommer* clone();

        virtual int nFaceCorners() { return 2; }
        virtual int* faceCorners() { int *fc =  new int[2]; fc[0] = nn[0]; fc[1] = nn[order-1]; return fc; }

        int* wetDofs(DofSetArray &, int *p=0);

        void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int pflag=0);
        void neumVectorDeriv(CoordSet&,ComplexVector&,
                                double,double,double,double, int, int pflag=0);
        void wetInterfaceVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int,int);

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
                                                        double *d);
        void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);
        FullSquareMatrixC sommer2Matrix(CoordSet&, complex<double> *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return exp(ComplexD(0.0,M_PI/4.0))/sqrt(8.0*M_PI*k); }

	void getNormal(CoordSet&, double [3]);

        void markDofs(DofSetArray &);

};
#endif
