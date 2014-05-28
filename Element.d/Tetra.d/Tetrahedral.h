#ifndef _TETRAHEDRAL_H_
#define _TETRAHEDRAL_H_

#include <Element.d/Element.h>
#include <Element.d/Tetra.d/TetraElementTemplate.hpp>
class TetCorotator;
class NLMaterial;

class Tetrahedral: public Element,
                   public TetraElementTemplate<double> {

	int nn[4];
        double *nodeTemperatures;
        double  *cCoefs;  // HB 4-15-05
        double  *cFrame;  // HB 4-15-05

        TetCorotator *tetCorotator;
        NLMaterial *mat;
public:
	Tetrahedral(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *mel, int cmflg=1);
        double weight(CoordSet& cs, double *gravityAcceleration);
	double           getMass(CoordSet& cs);
        void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);
        void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force,
                              int glflag, GeomState *gs=0);

        void             getVonMises
                                (Vector &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface=0,
                                 double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void getVonMisesNodalCoordinateSensitivity
                                (GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                 int senMethod=1, double* ndTemps=0, int avgnum=1, double ylayer=0, double zlayer=0);

        void             getAllStress
                                (FullM &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface=0,
                                 double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
        Corotator *getCorotator(CoordSet &cs, double *kel = 0, int d1 = 2, int d2 = 2);

	// DEC
	PrioInfo examine(int sub, MultiFront *);
        int nDecFaces() { return 3;}
        int getDecFace(int iFace, int *fn);

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame)
          { cCoefs = coefs; cFrame = frame; } // PJSA 3-30-05

        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Tetrahedral el.\n"); return (double *) 0;
        }
	
        //HB 05-26-05 
        void getVonMisesAniso(Vector &stress, Vector &weight,CoordSet &cs,
          	              Vector &elDisp,int strInd,int surface=0,
        	              double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);


        void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
		               Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

        void setMaterial(NLMaterial *);
        int numStates();
};
#endif
