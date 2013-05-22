// ---------------------------------------------------------------------
// HB - 01-15-06
// ---------------------------------------------------------------------
// 15 nodes wedge element
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#ifndef _PENTA15_H_
#define _PENTA15_H_

#include <Element.d/Element.h>
class Corotator;
class NLMaterial;

class Penta15: public Element {

	int nn[15];
        double *C; 

        Corotator* penta15Corotator; 

        double  *cCoefs;  
        double  *cFrame; 
        NLMaterial *mat;

public:
	Penta15(int*);

        Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

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
                 "              for Penta15 el.\n"); return (double *) 0;
        }

        void setMaterial(NLMaterial *);
        int numStates();
};
#endif

