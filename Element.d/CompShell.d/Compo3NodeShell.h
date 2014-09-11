#ifndef _COMPO3NODESHELL_H_
#define _COMPO3NODESHELL_H_

#include <Element.d/Element.h>

class GeomState;

class Compo3NodeShell : public Element {

        int      nn[3];
        int      type;
        int numLayers;
        int     (*idlay)[5];
        double  *layData;
        double  *cCoefs;
        double  *cFrame;
        PressureBCond *pbc;
        static bool Wzero_density, Wthermal_force;

public:
	Compo3NodeShell(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs,double *mel,int cmflg=1);

        void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                     GeomState *gs);
        void getGravityForceSensitivityWRTthickness(CoordSet& cs, double *gravityAcceleration,
                                                    Vector& gravityForceSensitivity, int gravflg, GeomState *geomState);

        void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface=0,
                         double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface=0,
                         double *ndTemps=0);

        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame);
        double * setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta);

	void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int  numDofs();
	int  getTopNumber();

        int  numNodes();
        int* nodes(int * = 0);
	Corotator *getCorotator(CoordSet &, double *,int,int);
        double getMass(CoordSet &);
        double getMassSensitivityWRTthickness(CoordSet &);
        double weight(CoordSet&, double *);
        double weightDerivativeWRTthickness(CoordSet&, double *, int=1);

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs=0);
        void getFlLoad(CoordSet &, const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs=0);

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);
        void getThermalForce(CoordSet&, Vector&, Vector&, int glflag, 
	                     GeomState *gs=0);

        virtual int getCompositeLayer() { return numLayers;  }

	virtual double * getMidPoint(CoordSet &);
        virtual double * getCompositeData(int nl) { return layData+(nl*9); }
        virtual double * getCompositeFrame()      { return cFrame;  }

#ifdef STRUCTOPT

protected:   
	double  *layGrad;
	
public:	
        int chkOptInf(CoordSet&);

        double getGradMass(CoordSet& cs, CoordSet& dcs);
	
	void setCompositeGrad( double *gradval ) {layGrad=gradval;}

        void getGradVonMises(Vector &dstress, Vector &weight, 
	                     CoordSet &cs, CoordSet &dcs, 
			     Vector &elDisp, Vector &elGrad,
			     int strInd, int surface);

        void gradMassMatrix(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dmel);

        void gradstiffness(CoordSet& cs, CoordSet& gradcs, 
                           FullSquareMatrix & gradkarray, int flg=1);
#endif
	
        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);
        int nDecFaces() { return 1; }
        int getDecFace(int iFace, int *fn) { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

        int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }

	bool hasRot() { return true; }

        int getMassType() { return 0; } // lumped only

};
#endif

