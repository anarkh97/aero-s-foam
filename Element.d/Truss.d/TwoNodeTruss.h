#ifndef _TWONODETRUSS_H_
#define _TWONODETRUSS_H_

#include <Element.d/Element.h>
#include <Element.d/Truss.d/TrussElementTemplate.hpp>

class TwoNodeTruss : public virtual Element,
                     public TrussElementTemplate<double> 
{
        int nn[2];
        double preload;
public:
	TwoNodeTruss(int*);
        Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
        int getMassType() { return 2; } // both consistent and lumped
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double getMass(CoordSet&);
        void getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx);
        void getLengthNodalCoordinateSensitivity(CoordSet &cs, Vector &dLengthdx);
        double weight(CoordSet& cs, double *gravityAcceleration);
        void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration);
        void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg, GeomState *gs);
        void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                          GenFullM<double> &dGfdx, int gravflg, GeomState *geomState=0);
        void getIntrnForce(Vector &elForce, CoordSet& cs,
 	                   double *elDisp, int forceIndex, double *ndTemps=0);
        void getVonMises(Vector& stress, Vector& weight,CoordSet &cs,
	                 Vector& elDisp, int strInd,int surface=0, 
	                 double *ndTemps=0,double ylayer=0.0, double zlayer=0.0, int avgnum=1);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, 
                                                Vector &elDisp, int strInd, int surface,
                                                double *ndTemps, int avgnum, double ylayer, double zlayer);
        void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                   double *ndTemps, int avgnum, double ylayer, double zlayer);
        void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int numDofs();
        int numNodes();
        int* nodes(int * = 0);
        Corotator* getCorotator(CoordSet &cs, double *kel,int,int);
        int getTopNumber();
        void getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &ThermalForce, int glflag, 
                             GeomState *gs);
        void setPreLoad(std::vector<double> &load);
        std::vector<double> getPreLoad();
	bool isSafe() { return false; }
	bool isStart() {return false; }
#ifndef SALINAS
        PrioInfo examine(int sub, MultiFront *);
#endif
private:
        void buildBarFrame(CoordSet&, double xg[3][3], double xl[3][3]);
};
#endif
