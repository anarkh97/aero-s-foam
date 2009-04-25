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

public:
	Compo3NodeShell(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs,double *mel,int cmflg=1);

        void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                     GeomState *gs);

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

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs=0);
        void getFlLoad(CoordSet &, const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs=0);
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs, int cflg);

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

	bool hasRot() { return true; }
	//double weight() { return 3; }
        //double trueWeight() { return 3; }
};
#endif

