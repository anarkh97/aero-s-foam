#ifndef _TENNODETETRAHEDRAL_H_
#define _TENNODETETRAHEDRAL_H_

#include <Element.d/Element.h>

class TenNodeTetrahedral: public Element {

	int nn[10];

        double  *cCoefs;  // HB 06-19-05
        double  *cFrame;  // HB 06-19-05
                                                                                                                                        
public:
	TenNodeTetrahedral(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);
        void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);
        void             getThermalForce(CoordSet &cs,Vector &nT,Vector &force, 
                                         int glflag, GeomState *gs=0);

        void             getVonMises
                                (Vector &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface=0,
                                 double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

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

	PrioInfo examine(int sub, MultiFront *);

	//double weight() { return 3.0; }
	//double trueWeight() { return 3.0; }
	//int facelist(PolygonSet &pgset, int *facelist=0);

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame)
          { cCoefs = coefs; cFrame = frame; } // HB 06-19-05 
                                                                                                                                        
        double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
        { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for Tetrahedral el.\n"); return (double *) 0;
        }

        //HB 06-11-05 
        void getVonMisesAniso(Vector &stress, Vector &weight,CoordSet &cs,
          	              Vector &elDisp,int strInd,int surface=0,
        	              double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                               Vector &elDisp, int strInd, int surface=0,
                               double *ndTemps=0);

};
#endif

