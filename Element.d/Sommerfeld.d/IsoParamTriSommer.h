#ifndef _ISOPARAMTRISOMMER_H_
#define _ISOPARAMTRISOMMER_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamTriSommer: public SommerElement {

	int *nn;
        int order;
public:
	IsoParamTriSommer(int, int*, Element *_el = 0, int etype = 11);

        int numNodes() {return (order*(order+1))/2;}
        int getNode(int nd) { return nn[nd]; }
	int* getNodes() { return nn; }
        int numDofs() { return (order*(order+1))/2; }
        int numSolidDofs() { return 3*(order*(order+1))/2; }
        int numWetDofs() { return 2*(order*(order+1)); }
        int dim() { return 3; }
        int* dofs(DofSetArray &, int *p=0);
        int* wetDofs(DofSetArray &, int *p=0);
        int* solidDofs(DofSetArray &, int *p=0);
        virtual IsoParamTriSommer* clone();

        virtual int nFaceCorners() { return 3; }
        virtual int* faceCorners() { int *fc =  new int[3]; fc[0] = nn[0];
          fc[1] = nn[order-1]; fc[2] = nn[(order*(order+1))/2-1];  return fc; }

        void flipNormal();

        void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double, int pflag = 0);
        void neumVectorDeriv(CoordSet& cs, ComplexVector& cv, double k,
                                    double dx, double dy, double dz, int n,
                                    int pflag=0);
        void wetInterfaceVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int,int);
        void wetInterfaceVector(CoordSet&,ComplexVector&,
                                complex<double> (*)[3],complex<double>*);

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
                                                        double *d);
        void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);
        virtual void ffp(CoordSet &cs, int numFFP, double *dirFFP,
                         complex<double> *sol, complex<double> *ffpv, bool direction);

        FullSquareMatrixC sommer2Matrix(CoordSet&, complex<double> *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return ComplexD(0.25/M_PI, 0.0); }

	void getNormal(CoordSet&, double [3]);

        void markDofs(DofSetArray &);
};
#endif
