#ifndef _THREENODESHELL_H_
#define _THREENODESHELL_H_

#include  <Element.d/Element.h>
#include	<Element.d/Shell.d/ShellElementSemiTemplate.hpp>
#include        <Driver.d/MultiFront.h>
class GeomState;
class Shell3Corotator;

class ThreeNodeShell : virtual public Element,
                       public ShellElementSemiTemplate<double> {
protected:
	int nn[3];
	double w;
        Shell3Corotator *corot;
        PressureBCond *pbc;
public:
	ThreeNodeShell(int*, double _w=3);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);
        double weight(CoordSet&, double *, int);
        double weightDerivativeWRTthickness(CoordSet&, double *, int, int=1);
        void   getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
                                         GeomState *gs);
        void   getVonMises(Vector& stress, Vector& weight,
                           CoordSet &cs, Vector& elDisp, int strInd,int surface=0,
                           double *ndTemps=0,double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                             int senMethod, double *, int avgnum, double ylayer, double zlayer);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                int senMethod, double *ndTemps, int avgnum, double ylayer, double zlayer);
        void             getAllStress(FullM& stress, Vector& weight,
                         CoordSet &cs, Vector& elDisp, int strInd,int surface=0,
                                 double *ndTemps=0);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        Corotator *getCorotator(CoordSet &, double *, int , int);

        void             computeDisp(CoordSet&, State &, const InterpPoint &,
                                     double*, GeomState *gs);
        void printInfo(CoordSet&, State &, double[2]);
	void		 getFlLoad(CoordSet &, const InterpPoint &,
                                   double *flF, double *resF, GeomState *gs=0);

	int getTopNumber();

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0 );

	void getThermalForce(CoordSet& cs, Vector& ndTemps,Vector &elThermalForce, 
	                     int glfag, GeomState *gs=0);

        bool isShell() { return true; }

        int getMassType() { return 0; } // lumped only

#ifdef STRUCTOPT

        int chkOptInf(CoordSet&);

        double getGradMass(CoordSet& cs, CoordSet& dcs);
	
        void getGradVonMises(Vector &dstress, Vector &weight, 
	                     CoordSet &cs, CoordSet &dcs, 
			     Vector &elDisp, Vector &elGrad,
			     int strInd, int surface);

        void gradMassMatrix(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dmel);

        void gradstiffness(CoordSet& cs, CoordSet& gradcs, 
                           FullSquareMatrix & gradkarray, int flg=1);

#endif
	// DEC
	bool hasRot() {return true;}
	PrioInfo examine(int sub, MultiFront *mf);

};
#endif

