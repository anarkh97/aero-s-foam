#ifndef _QUADSOMMERBC_H_
#define _QUADSOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

class QuadSommerBC: public SommerElement {

	int nn[4];
public:
	QuadSommerBC(int, int, int, int, Element *_el = 0, int eType = 4);

        int numNodes() {return 4;}
        int getNode(int nd) { return nn[nd]; }
        int* getNodes() { return nn; }
        int  numDofs();
        int dim() { return 3; }
        int* dofs(DofSetArray &, int *p=0);
        virtual QuadSommerBC* clone();

        FullSquareMatrix sommerMatrix(CoordSet&, double *);
        FullSquareMatrix turkelMatrix(CoordSet&, double *);
        FullSquareMatrix refinedSommerMatrix(CoordSet&, double *);
        //FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
        FullSquareMatrix HSommerMatrix(CoordSet&, double *);
        //FullSquareMatrix HKSommerMatrix(CoordSet&, double *);

        ComplexD ffpCoef(double k) { return ComplexD(0.25/M_PI, 0.0); }

	void getNormal(CoordSet&, double[3]);
        void markDofs(DofSetArray &);

private:
	void getLocalCoordinates (CoordSet&, double xx[4], double yy[4], double zz[4]);
        void SurfaceRefinement(int nNo, double* x, double* y, double* z, double* xx, double* yy, double* zz);
        void GaussCoordinates(int Ngp, double* Pg, double* weight);
};
#endif

