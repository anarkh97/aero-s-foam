#ifndef _TRIANGLESOMMERBC_H_
#define _TRIANGLESOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class TriangleSommerBC: public SommerElement {

	int nn[3];
public:
	TriangleSommerBC(int, int, int, Element *_el = 0, int eType = 3);

        int numNodes() {return 3;}
        int getNode(int nd) { return nn[nd]; }
        int* getNodes() { return nn; }
        int numDofs(); 
        int dim() { return 3; }
        int* dofs(DofSetArray &, int *p=0);
        void markDofs(DofSetArray &);
        virtual TriangleSommerBC* clone();

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);
        FullSquareMatrix refinedSommerMatrix(CoordSet&, double *);
        //FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
        FullSquareMatrix HSommerMatrix(CoordSet&, double *);
        //FullSquareMatrix HKSommerMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return ComplexD(0.25/M_PI, 0.0); }

        void getNormal(CoordSet&, double[3]);
        double getSize(CoordSet&);

        void BT2(CoordSet& cs, double *e, double *f, double *g,
                 double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d);
        void BT2n(CoordSet& cs, double *e, double *f, double *g,
                 double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n);
        void sphereBT2(CoordSet& cs, double r, double k, ComplexD *d);


private:
        void get_basis(int, int, int, double (*)[3], double*, double*); 
        double getArea(CoordSet&, int*); 
	void getLocalCoordinates (CoordSet&, double xx[3], double yy[3], double
zz[3]);
        void SurfaceRefinement(int nNo, double* x, double* y, double* z, double* xx, double* yy, double* zz);
        void GaussCoordinates(int Ngp, double* uPg, double* vPg, double* weight);

};
#endif

