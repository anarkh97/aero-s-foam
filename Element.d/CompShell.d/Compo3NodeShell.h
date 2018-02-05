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

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
        FullSquareMatrix massMatrix(const CoordSet& cs,double *mel,int cmflg=1) const;

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
        void getCFrame(CoordSet &cs, double cFrame[3][3]) const;

        void markDofs(DofSetArray &) override;
        int* dofs(DofSetArray &, int *p) override;
         int numDofs() const override;
        int getTopNumber() override;

        int numNodes() const;
        int* nodes(int * = 0) const override;
        Corotator *getCorotator(CoordSet &, double *,int,int);
        double getMass(const CoordSet &) const;
        double getMassThicknessSensitivity(CoordSet &);

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

        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *) override;
        int nDecFaces() const override { return 1; }
        int getDecFace(int iFace, int *fn) { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

        int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }

	bool hasRot() { return true; }

        int getMassType() const override { return 0; } // lumped only

};
#endif

