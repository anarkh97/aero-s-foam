#ifndef _SHEARPANEL_H_
#define _SHEARPANEL_H_

#include <Element.d/Element.h>
#include <Element.d/Shear.d/ShearPanelTemplate.hpp>

class ShearPanel: public Element,
                  public ShearPanelTemplate<double> {

	int nn[4];
public:
	ShearPanel(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double getMass(CoordSet&);
        double getMassSensitivityWRTthickness(CoordSet&);
        double weight(CoordSet&, double *);
        double weightDerivativeWRTthickness(CoordSet&, double *, int=1);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *geomState);
        void getGravityForceSensitivityWRTthickness(CoordSet&,double *gravity, int senMethod, Vector& f, int gravflg,
	                                                  GeomState *geomState);

        void             getVonMises (Vector &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0, 
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                int senMethod, double *ndTemps, int avgnum, double ylayer, double zlayer);

        void             getAllStress(FullM &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	PrioInfo examine(int sub, MultiFront *);
	int getTopNumber();
	bool hasRot() {return true;}

        int getMassType() { return 0; } // lumped only
};
#endif

