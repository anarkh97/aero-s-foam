// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 32 nodes brick element
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#ifndef _BRICK32_H_
#define _BRICK32_H_

#include <Element.d/Element.h>

class Brick32: public Element {

	int nn[32];
        double *C; 

        Corotator* brick32Corotator; 

        double  *cCoefs;  
        double  *cFrame; 

public:
	Brick32(int*);

        Element *clone();

	void renum(int *);

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void		 getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);

        void             getVonMises(Vector &stress, Vector &weight, 
                                     CoordSet &cs, Vector &elDisp, 
                                     int strInd, int surface=0,
                                     double *ndTemps=0,
				     double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void             getAllStress(FullM &stress, Vector &weight,
                                  CoordSet &cs, Vector &elDisp,
                                  int strInd, int surface=0,
                                  double *ndTemps=0);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
        int numTopNodes();
        Corotator *getCorotator(CoordSet &cs, double *kel, int , int );

       void getThermalForce(CoordSet &cs, Vector &ndTemps,
                            Vector &elementThermalForce, int glflag, 
			    GeomState *geomState);

       PrioInfo examine(int sub, MultiFront *); // dec
       //double weight() { return 3; } // should be 12 for brick20
       //double trueWeight() { return 3; } // should be 12 for brick20

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame) 
          { cCoefs = coefs; cFrame = frame; } 

        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Hexahedral el.\n"); return (double *) 0;
        }

};
#endif

