#ifndef _MEMBRANE_H_
#define _MEMBRANE_H_

#include <Element.d/Element.h>

class Membrane : public Element
{

        int nn[3];
public:
        Membrane(int*);

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
        int getTopNumber();
        Corotator* getCorotator(CoordSet &cs, double *kel, int fitAlg, int);

        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);

        // Miscellaneous
        bool hasRot() { return true; }
        int  getMassType() { return 0; } // lumped only

        // NEW STRUCTOPT 
        double getMassThicknessSensitivity(CoordSet& cs);
        void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet&, double *gravityAcceleration,
                                                int senMethod = 1);
        void getGravityForceThicknessSensitivity(CoordSet&, double *gravity, int senMethod, Vector&, int gravflg,
                                                 GeomState *gs = 0);
        void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration, int senMethod,
                                                       GenFullM<double> &dGfdx, int gravflg, GeomState *gs = 0);
        void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg = 1,
                                              int senMethod = 0);
        void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs, int senMethod = 0);
        void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                   Vector &elDisp, int strInd, int surface, int senMethod = 1,
                                                   double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
        void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp,
                                             int strInd, int surface, int senMethod = 1, double *ndTemps = 0, int avgnum = 1,
                                             double ylayer = 0, double zlayer = 0);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs,
                                                Vector &elDisp, int strInd, int surface, int senMethod = 1,
                                                double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
};
#endif
