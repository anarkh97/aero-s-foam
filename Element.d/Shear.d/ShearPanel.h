#ifndef _SHEARPANEL_H_
#define _SHEARPANEL_H_

#include <Element.d/Element.h>
#include <Element.d/Shear.d/ShearPanelTemplate.hpp>

class ShearPanel: public Element,
                  public ShearPanelTemplate<double> {

	int nn[4];
public:
	ShearPanel(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double getMass(CoordSet&);
        double getMassThicknessSensitivity(CoordSet&);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *geomState);
        void getGravityForceThicknessSensitivity(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                               GeomState *geomState);

        void             getVonMises (Vector &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0, 
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                                CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                double *ndTemps, int avgnum, double ylayer, double zlayer);

        void             getAllStress(FullM &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;
	bool hasRot() {return true;}

        int getMassType() { return 0; } // lumped only
};
#endif

