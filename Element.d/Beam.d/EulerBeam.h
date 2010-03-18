#ifndef _EULERBEAM_H_
#define _EULERBEAM_H_

#include <Element.d/Element.h>

class GeomState;

class EulerBeam : public Element 
{
        EFrame *elemframe;
        int nn[3];
	double *offset;
        double c0[3][3];
public:

	EulerBeam(int*);

	Element *clone();

	void renum(int *);

	void setFrame(EFrame *ef) { elemframe = ef; }
        void buildFrame(CoordSet&);

	FullSquareMatrix stiffness(CoordSet&, double *kel,int flg=1);
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);
	void             getGravityForce(CoordSet&, double *g, Vector &f, int gravflg, 
	                                 GeomState *gs);
        void             getIntrnForce(Vector& elForce, CoordSet& cs,
				       double* elDisp,int forceIndex, double *ndTemps=0);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

        Corotator *getCorotator(CoordSet &, double*, int, int);
	int getTopNumber();
	void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs, int cflg);

        void getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &ThermalForce, int glflag, GeomState *gs);
        void computeDisp(CoordSet &cs, State &state, const InterpPoint &,
                         double *res, GeomState *gs);
        void getFlLoad(CoordSet &cs,  const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs);
	void setOffset(double *o) { offset = o; }
	
	void getVonMises(Vector& stress, Vector& weight,CoordSet &cs, Vector& elDisp, 
	                 int strInd,int surface=0, double *ndTemps=0, 
			 double ylayer=0.0, double zlayer=0.0, int avgnum=0);

	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *);
	bool hasRot() { return true; }
   
        int getMassType() { return 2; } // both consistent and lumped

private:

       void    getLength(CoordSet&, double &length);
       void    updTransMatrix(CoordSet&, GeomState *gs, double t[3][3], double &len);
       void    offsetAxis(FullSquareMatrix& mat);
       
};
#endif

