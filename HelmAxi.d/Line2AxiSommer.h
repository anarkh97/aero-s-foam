#ifndef _LINE2AXISOMMER_H_
#define _LINE2AXISOMMER_H_ 

#include <Utils.d/MyComplex.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>

template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

class Line2AxiSommer: public SommerElement {

	int nn[3];
	int type;
        double surfR0, surfZ0;

public:
	Line2AxiSommer(int, int, int);

        void setType(int t);
        void setSurf(double aR, double aZ);

        void renum(int *);

        FullSquareMatrix sommerMatrix(CoordSet&);
        FullSquareMatrix sommerMatrix(CoordSet&, double *);

        FullSquareMatrixC turkelMatrix(CoordSet&, double, int);
        FullSquareMatrixC turkelMatrix(CoordSet&, double, int, DComplex *);

        FullSquareMatrix interfMatrixConsistent(CoordSet&);
        FullSquareMatrix interfMatrixConsistent(CoordSet&, double*);
        FullSquareMatrix interfMatrixLumped(CoordSet&);
        FullSquareMatrix interfMatrixLumped(CoordSet&, double*);

        int* dofs(DofSetArray &, int *p=0);
        int  numDofs();
	int* nodes(int* = 0);
        int  numNodes();

        void ffpAxiNeum(int, DComplex *, CoordSet &, DComplex **, double,
                     double(*)[3], double*, int);
        void ffpAxiDir(int, DComplex *, CoordSet &, DComplex **, double,
                     double(*)[3], double*, int);

};

#endif
