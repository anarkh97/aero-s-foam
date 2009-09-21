#ifndef	_TIMOSHENKOBEAM_H_
#define	_TIMOSHENKOBEAM_H_

#include        <Element.d/Element.h>

class GeomState;

class TimoshenkoBeam : public Element {
        EFrame *elemframe;

        int nn[3];

  	double* iniOr;	
	
public:
        TimoshenkoBeam(int*);

	Element *clone();

	void renum(int *);

        void setFrame(EFrame *ef) { elemframe = ef; }
        int buildFrame(CoordSet&);

	FullSquareMatrix stiffness(CoordSet& cs, double *kel, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

	double           getMass(CoordSet& cs);
	void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);
        void             getIntrnForce(Vector& elForce,CoordSet& cs,
				       double *elDisp,int forceIndex, double *ndTemps);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        Corotator *getCorotator(CoordSet &, double *, int,int);
	int getTopNumber();
	void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs, int cflg);
	
	void getThermalForce(CoordSet &, Vector &, Vector &, int, GeomState *geomState=0);
     	void getVonMises(Vector& stress, Vector& weight,CoordSet &cs, Vector& elDisp, 
	                 int strInd,int surface=0, double *ndTemps=0,
			 double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        int getMassType() { return 0; } // lumped only
private:
       
        void updTransMatrix(CoordSet&, GeomState *gs, double t[3][3], double &len);

#ifdef STRUCTOPT
    protected:
       EFrame *delemframe;
    
    public:	

        int chkOptInf(CoordSet&);

        double getGradMass(CoordSet&, CoordSet&);

        double gradFrame(CoordSet&, CoordSet&);

        void setdFrame(EFrame *def) { delemframe = def; }

        void getGradIntrnForce(Vector&, CoordSet&, CoordSet&,
	                       double*, double*, int , double*);

        void getGradGravityForce(CoordSet&, CoordSet&, double*, Vector&);
        
	void gradMassMatrix(CoordSet &, CoordSet &, FullSquareMatrix &); 
	
	void gradstiffness (CoordSet&, CoordSet&, FullSquareMatrix &, int=1);
	
        void updFrame(CoordSet&);

#endif

	// Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);
        // isStart indicates if an element is suitable to
        // be a germination center for a subdomain (bars are not)
        bool isStart() { return true; }  
	bool hasRot() { return true; }

};
#endif
