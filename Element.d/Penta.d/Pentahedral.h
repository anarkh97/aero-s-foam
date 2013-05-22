#ifndef _PENTAHEDRAL_H_
#define _PENTAHEDRAL_H_

#include <Element.d/Element.h>

class PentaCorotator;
class NLMaterial;

class Pentahedral: public Element {

	int nn[6];
        double  *cCoefs;  // HB 4-15-05
        double  *cFrame;  // HB 4-15-05
                                                                                                                  
        PentaCorotator* pentaCorotator;
        NLMaterial *mat;
public:
	Pentahedral(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

	FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);
        void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *geomState);

        void getVonMises(Vector &stress, Vector &weight,CoordSet &cs,
			 Vector &elDisp,int strInd,int surface=0,
                                 double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void getAllStress(FullM &stress, Vector &weight,CoordSet &cs,
			 Vector &elDisp,int strInd,int surface=0,
                                 double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();

        // HB (09-21-03): implement thermal forces. MUST BE VALIDATED !!!
       void getThermalForce(CoordSet &cs, Vector &ndTemps,
                            Vector &elementThermalForce, int glflag, GeomState *geomState);

       PrioInfo examine(int sub, MultiFront *);
       //double weight() { return 2; } // should be 12 for penta 15
       //double trueWeight() { return 2; } // should be 12 for penta 15

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame)
          { cCoefs = coefs; cFrame = frame; } // PJSA 3-30-05

        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Pentahedral el.\n"); return (double *) 0;
        }

       //HB: add "corotator" methods
       Corotator *getCorotator(CoordSet &cs, double *kel, int , int );
       
      //HB 05-26-0 
      void getVonMisesAniso(Vector &stress, Vector &weight,CoordSet &cs,
        	            Vector &elDisp,int strInd,int surface=0,
        	            double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);


     void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
		            Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

     void setMaterial(NLMaterial *);
     int numStates();
};
#endif

