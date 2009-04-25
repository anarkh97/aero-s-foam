#ifndef _TRIANGLE6SOMMERBC_H_
#define _TRIANGLE6SOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

typedef double Matrix22[2][2];

class Triangle6SommerBC: public SommerElement {

	int nn[6];
public:
	Triangle6SommerBC(int, int, int, int, int, int, Element *_el = 0, int etype = 6);

        int numNodes() {return 6;}
        int getNode(int nd) { return nn[nd]; }
        int* getNodes() { return nn; }
        int  numDofs();
        int dim() { return 3; }
        int* dofs(DofSetArray &, int *p=0);
        void flipNormal();
        virtual Triangle6SommerBC* clone();

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        void sommerMatrixEllipsoid(CoordSet &cs, double kappa,  double H[3], double K[3], ComplexD *d);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return ComplexD(0.25/M_PI, 0.0); }

	void getNormal(CoordSet&, double[3]);

        void BT2(CoordSet& cs, double *e, double *f, double *g,
                 double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d);
        void sphereBT2(CoordSet& cs, double r, double k, ComplexD *d);
        void ellipsoidBT2(CoordSet& cs, double a, double b, double k, ComplexD *d);

//        void renum (int *);
        void markDofs(DofSetArray &);
//        FullSquareMatrix  stiffness(CoordSet&, double *d, int flg = 1);
//        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
//        int* nodes(int * = 0);

//        bool isSommerElement() { return true; }
private:
        static double tri6_derivatives[12][6][2];
        static double tri6_values[12][6];
        static double tri6_weights[12];
        static double tri6_coord[12][2];
        static double tri3_values[12][3];
        void computedxdxi(double x[6], double y[6], Matrix22 *dxdxi, double *det);
        void get_basis(int, int, int, double (*)[3], double*, double*);
        double getArea(CoordSet&, int*);
	void getLocalCoordinates (CoordSet&, double xx[6], double yy[6], double zz[6]);
        void getLocalCoordinatesNew (CoordSet&, double T[3][3], double xi1[3], double xi2[3], double xi3[3]);


};
#endif
