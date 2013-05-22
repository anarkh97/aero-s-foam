#ifndef _TWONODETRUSSF_H_
#define _TWONODETRUSSF_H_

#include <Element.d/Element.h>
class BarFCorotator;

class TwoNodeTrussF : public virtual Element {

        int nn[2];
        double preload;
        BarFCorotator *myCorot;
public:
	TwoNodeTrussF(int*);
        Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double getMass(CoordSet&);
        void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
	                     GeomState *gs);
        void getIntrnForce(Vector &elForce, CoordSet& cs,
 	                   double *elDisp, int forceIndex, double *ndTemps=0);
        void getVonMises(Vector& stress, Vector& weight,CoordSet &cs,
	                 Vector& elDisp, int strInd,int surface=0, 
	                 double *ndTemps=0,double ylayer=0.0, double zlayer=0.0, int avgnum=0);
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

private:

        void buildBarFrame(CoordSet&, double xg[3][3], double xl[3][3]);

#ifndef SALINAS
	PrioInfo examine(int sub, MultiFront *);
#endif
};
#endif
