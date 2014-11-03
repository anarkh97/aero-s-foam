#ifndef _THREENODESHELL_H_
#define _THREENODESHELL_H_

#include <Element.d/Element.h>
#include <Driver.d/MultiFront.h>

class GeomState;
class Shell3Corotator;

class ThreeNodeShell : public Element
{
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
        double getMass(CoordSet& cs);
        void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg,
                             GeomState *gs);
        void getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
                         Vector& elDisp, int strInd, int surface=0,
                         double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void getAllStress(FullM& stress, Vector& weight, CoordSet& cs,
                          Vector& elDisp, int strInd, int surface=0,
                          double *ndTemps=0);

        void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int numDofs();

        int numNodes();
        int* nodes(int * = 0);
        Corotator *getCorotator(CoordSet &, double *, int , int);

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs);
        void printInfo(CoordSet&, State &, double[2]);
        void getFlLoad(CoordSet &, const InterpPoint &,
                       double *flF, double *resF, GeomState *gs=0);

        int getTopNumber();
        int nDecFaces() { return 1;}
        int getDecFace(int iFace, int *fn) { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

        int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0 );

        void getThermalForce(CoordSet& cs, Vector& ndTemps,Vector &elThermalForce, 
                             int glfag, GeomState *gs=0);

        bool isShell() { return true; }

        int getMassType() { return 0; } // lumped only

        // DEC
        bool hasRot() {return true;}
        PrioInfo examine(int sub, MultiFront *mf);

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
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs,
                                                Vector &elDisp, int strInd, int surface, 
                                                double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

};
#endif

