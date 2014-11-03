#ifndef	_TIMOSHENKOBEAM_H_
#define	_TIMOSHENKOBEAM_H_

#include <Element.d/Element.h>
#include <Element.d/Beam.d/BeamElementTemplate.hpp>

class GeomState;

class TimoshenkoBeam : public Element,
                       public BeamElementTemplate<double> {
        EFrame *elemframe;
        double oeframe[3][3];
        bool myElemFrame;
        int nn[3];
        double* iniOr;	
        PressureBCond *pbc;
	
public:
        explicit TimoshenkoBeam(int*);
        virtual ~TimoshenkoBeam();

        TimoshenkoBeam *clone();

        void renum(int *);
        void renum(EleRenumMap&);

        void setFrame(EFrame *ef) { elemframe = ef; myElemFrame = false; }
        const EFrame *getFrame() const { return elemframe; }
        void buildFrame(CoordSet&);

        FullSquareMatrix stiffness(CoordSet& cs, double *kel, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

        void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
        double getMass(CoordSet& cs);
        void getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx);
        double weight(CoordSet& cs, double *gravityAcceleration);
        void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration);
        void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg, GeomState *gs);
        void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                       GenFullM<double> &dGfdx, int gravflg, GeomState *geomState=0);
        void getIntrnForce(Vector& elForce,CoordSet& cs,
                           double *elDisp,int forceIndex, double *ndTemps);

        void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int numDofs();

        int numNodes();
        int* nodes(int * = 0);
        Corotator *getCorotator(CoordSet &, double *, int,int);
        int getTopNumber();

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);

        void getThermalForce(CoordSet &, Vector &, Vector &, int, GeomState *geomState=0);
        void getVonMises(Vector& stress, Vector& weight,CoordSet &cs, Vector& elDisp, 
                         int strInd,int surface=0, double *ndTemps=0,
                         double ylayer=0.0, double zlayer=0.0, int avgnum=1);
#ifdef USE_EIGEN3
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp,
                                                int strInd, int surface, double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

        void getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight, CoordSet &cs, ComplexVector &elDisp,
                                                int strInd, int surface, double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
        void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                   double *ndTemps=0, int avgnum=1, double ylayer=0, double zlayer=0);
#endif

        int getMassType() { return 0; } // lumped only
private:
        TimoshenkoBeam(const TimoshenkoBeam &);
        TimoshenkoBeam &operator=(const TimoshenkoBeam &);

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
	      bool hasRot() { return true; }
};

#endif
