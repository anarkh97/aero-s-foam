#ifndef _ISOPARAMQUADSOMMER_H_
#define _ISOPARAMQUADSOMMER_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamQuadSommer: public SommerElement {

	int *nn;
        int order;
public:
	IsoParamQuadSommer(int, int*, Element *_el = 0, int eType = 10);

        int numNodes() {return order*order;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return order*order; }
        int numSolidDofs() { return 3*order*order; }
        int numWetDofs() { return 4*order*order; }
        int dim() { return 3; }
        int* dofs(DofSetArray &, int *p=0);
        virtual IsoParamQuadSommer* clone();

        virtual int nFaceCorners() { return 4; }
        virtual int* faceCorners() { int *fc =  new int[4]; fc[0] = nn[0]; fc[1] = nn[order-1]; fc[2] = nn[order*order-1]; fc[3] = nn[order*(order-1)];  return fc; }


        int* wetDofs(DofSetArray &, int *p=0);
        int* solidDofs(DofSetArray &, int *p=0);

        void flipNormal();

        void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double, int pflag=0);
        void neumVectorDeriv(CoordSet& cs, ComplexVector& cv, double k,
                                    double dx, double dy, double dz, int n,
                                    int pflag=0);
        void wetInterfaceVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int,int);

        void wetInterfaceVector(CoordSet&,ComplexVector&,
                                complex<double> (*)[3],complex<double>*);
        void wetInterfaceVectorDeriv(CoordSet&,ComplexVector&,
                                complex<double> (*)[3],complex<double>*,
                                complex<double>*, int);

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
                                                        double *d);
        void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);
        virtual void ffp(CoordSet &cs, int numFFP, double *dirFFP,
                         complex<double> *sol, complex<double> *ffpv);

        FullSquareMatrixC sommer2Matrix(CoordSet&, complex<double> *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return ComplexD(0.25/M_PI, 0.0); }

	void getNormal(CoordSet&, double [3]);

        void markDofs(DofSetArray &);
};
#endif
