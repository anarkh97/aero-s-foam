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
        Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
        FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
        double getMass(const CoordSet&) const;
        void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
	                     GeomState *gs);
        void getIntrnForce(Vector &elForce, CoordSet& cs,
 	                   double *elDisp, int forceIndex, double *ndTemps=0);
        void getVonMises(Vector& stress, Vector& weight,CoordSet &cs,
	                 Vector& elDisp, int strInd,int surface=0, 
	                 double *ndTemps=0,double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p) const override;
         int numDofs() const override;
        int numNodes() const override;
        int* nodes(int * = 0) const override;
        Corotator* getCorotator(CoordSet &cs, double *kel,int,int);
	int getTopNumber() override;
        void getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &ThermalForce, int glflag, 
                             GeomState *gs);
        void setPreLoad(std::vector<double> &load);
        std::vector<double> getPreLoad();
	bool isSafe() const override { return false; }
	bool isStart() {return false; }

private:

        void buildBarFrame(CoordSet&, double xg[3][3], double xl[3][3]);

#ifndef SALINAS
	PrioInfo examine(int sub, MultiFront *) override;
#endif
};
#endif
