#ifndef _TWENTY_BRICK_H_
#define _TWENTY_BRICK_H_

#include <Element.d/Element.h>
class Corotator;
class NLMaterial;

class Brick20: public Element {

	int nn[20];
        double *C; 
        double  *cCoefs;  // HB 06-19-05
        double  *cFrame;  // HB 06-19-05
        NLMaterial *mat;                                                                                                                             
public:
	Brick20(int*);

        Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void		 getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);

        void             getVonMises (Vector &stress, Vector &weight, 
                                  CoordSet &cs, Vector &elDisp, 
                                  int strInd, int surface=0,
                                  double *ndTemps=0,
				  double ylayer=0.0, double zlayer=0.0, int avgnum=0);


        void             getAllStress(FullM &stress, Vector &weight, 
                                  CoordSet &cs, Vector &elDisp, 
                                  int strInd, int surface=0,
                                  double *ndTemps=0);
        void             getVonMisesInt(CoordSet &, Vector &, double &,
                                        double &, int, double &, double &);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();
        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
	// HB (09-21-03): implement WARNING message
	void getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &elementThermalForce, int glflag, GeomState *geomState);

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame)
          { cCoefs = coefs; cFrame = frame; } // HB 06-19-05 
                                                                                                                             
        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Hexa20 el.\n"); return (double *) 0;
        }

        void setMaterial(NLMaterial *);
        int numStates();
        Corotator *getCorotator(CoordSet &cs, double *kel, int, int);
};
#endif

