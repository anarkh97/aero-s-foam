#ifndef _LINEAXISOMMER_H_
#define _LINEAXISOMMER_H_ 

#include <Utils.d/MyComplex.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>

template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

class LineAxiSommer: public SommerElement {

        int nn[2];
        int type;
        double surfR0, surfZ0;

public:

	LineAxiSommer(int, int);
   
        void setType(int t);
        void setSurf(double aR, double aZ);

        void renum(int *);
        void renum(EleRenumMap&);

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
