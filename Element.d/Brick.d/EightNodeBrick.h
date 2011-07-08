#ifndef _EIGHTNODEBRICK_H_
#define _EIGHTNODEBRICK_H_

#include <Element.d/Element.h>

class BrickCorotator;

class EightNodeBrick: virtual public Element {
 protected:
	int nn[8];
        double *C; 
        BrickCorotator *brickCorotator;
        double  *cCoefs;  // PJSA 3-20-05
        double  *cFrame;  // PJSA 3-20-05


public:
	EightNodeBrick(int*);

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

        void             getIntrnForce(Vector &elForce, CoordSet& cs,
                                       double *elDisp, int forceIndex, double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
        Corotator *getCorotator(CoordSet &cs, double *kel, int , int );

       // HB (09-21-03): implement WARNING message: thermal force NOT implemented !!!
       void getThermalForce(CoordSet &cs, Vector &ndTemps,
                            Vector &elementThermalForce, int glflag, 
			    GeomState *geomState);

       PrioInfo examine(int sub, MultiFront *); // dec
       //double weight() { return 3; } // should be 12 for brick20
       //double trueWeight() { return 3; } // should be 12 for brick20

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame) 
          { cCoefs = coefs; cFrame = frame; } // PJSA 3-30-05

        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Hexahedral el.\n"); return (double *) 0;
        }

        //HB 05-26-0 
        void getVonMisesAniso(Vector &stress, Vector &weight,CoordSet &cs,
          	              Vector &elDisp,int strInd,int surface=0,
        	              double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);


        void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
		               Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);


};
#endif

