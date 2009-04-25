#ifndef _LINESOMMERBC_H_
#define _LINESOMMERBC_H_ 

#include <Element.d/Sommerfeld.d/SommerElement.h>

#ifdef WINDOWS
#include <math.h>
#endif

class LineSommerBC: public SommerElement {

	 int nn[2];
public:
	 LineSommerBC(int, int, Element *_el = 0, int etype = 1);

         int numNodes() {return 2;}
         int getNode(int nd) { return nn[nd]; }
	 int* getNodes() { return nn; }
         int* dofs(DofSetArray &, int *p=0);
         int numDofs() { return 2; }
         int dim() { return 2; }
         virtual LineSommerBC* clone();

         FullSquareMatrix sommerMatrix(CoordSet&, double *);
         FullSquareMatrix turkelMatrix(CoordSet&, double *);

         void sommerVector(CoordSet&, ComplexVector&, ComplexVector&);
         void btVector(CoordSet&, ComplexVector&, ComplexVector&);

         void ffpDir(int, ComplexD*, CoordSet&,ComplexD*,ComplexD*,
                     double,double(*)[3],double*);
         ComplexD ffpCoef(double k) { 
           return exp(complex<double>(0.0,M_PI/4.0))/sqrt(8.0*M_PI*k); }

	 void getNormal(CoordSet&, double [3]);
	 double getSize(CoordSet&);
};
#endif
