#ifndef _FELIPPASHELL_H_
#define _FELIPPASHELL_H_

#ifdef USE_EIGEN3
#include <Element.d/Element.h>
#include <Corotational.d/Shell3Corotator.h>

template <typename doublereal> class ShellMaterial;

class FelippaShell : public Element, 
                     public Shell3Corotator
{
        static int sflg, tflg; // see comments in ShellElementTemplate.hpp
        int      nn[3];
        int      type;
        double  *cFrame;
        PressureBCond *pbc;
        ShellMaterial<double> *gpmat; // gauss points' material
        ShellMaterial<double> *nmat;  // nodes material

public:
	FelippaShell(int*);
        ~FelippaShell();

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);

        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);

        void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg,
	                     GeomState *gs);

        void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface = 0,
                         double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
                         int avgnum = 0);

        void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface = 0,
                         double *ndTemps = 0);

        void setProp(StructProp *p, bool myProp = false);
        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame);
        double * setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta);
        void getCFrame(CoordSet &cs, double cFrame[3][3]) const;
        void setMaterial(NLMaterial *);
        int numStates();

	void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int  numDofs();
	int  getTopNumber();

        int  numNodes();
        int* nodes(int * = 0);
        double getMass(CoordSet &);

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs = 0);
        void getFlLoad(CoordSet &, const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs = 0);

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);
        void getThermalForce(CoordSet&, Vector&, Vector &, int, GeomState * = 0);

        // Nonlinear
        Corotator* getCorotator(CoordSet&, double*, int, int);
        void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                              FullSquareMatrix &elK, double *f, double dt, double t);
        void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                              FullSquareMatrix &elK, double *f, double dt, double t);
        void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt = 0);
        bool checkElementDeletion(GeomState &);
        void initStates(double *);
        double getDissipatedEnergy(GeomState &, CoordSet &);
        void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                                 int &nlflag);
        void getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                           GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
                           double ylayer = 0, double zlayer = 0, int avgnum = 0, int measure = -1);
        void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
                            GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
                            int measure = -1);

        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);
        int nDecFaces() { return 1; }
        int getDecFace(int iFace, int *fn) { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

        // Miscellaneous
        int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }
	bool hasRot() { return true; }
        int getMassType() { return 0; } // lumped only

        // NEW STRUCTOPT 
        double getMassThicknessSensitivity(CoordSet& cs);
        void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet&, double *gravityAcceleration);
        void getGravityForceThicknessSensitivity(CoordSet&, double *gravity, Vector&, int gravflg,
                                                 GeomState *gs = 0);
        void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                       GenFullM<double> &dGfdx, int gravflg, GeomState *gs = 0);
        void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg = 1);
        void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
        void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                   Vector &elDisp, int strInd, int surface,
                                                   double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
        void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp,
                                             int strInd, int surface, double *ndTemps = 0, int avgnum = 1,
                                             double ylayer = 0, double zlayer = 0);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet &cs,
                                                Vector &elDisp, int strInd, int surface,
                                                double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
        void getInternalForceThicknessSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs, Vector &dFintdThick,
                                                  double dt, double t);
  private:
        void getVonMisesImpl(Vector &stress, Vector &weight, CoordSet &cs,
                             Vector &elDisp, int strInd, int surface,
                             double *ndTemps, double ylayer, double zlayer,
                             int flag, double *staten = 0, double *statenp = 0);
        void getAllStressImpl(FullM &stress, Vector &weight, CoordSet &cs,
                              Vector &elDisp, int strInd, int surface,
                              double *ndTemps, int flag, double *staten = 0,
                              double *statenp = 0);
};
#endif
#endif

