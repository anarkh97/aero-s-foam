#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>
#include <Element.d/Quad4.d/QuadElementTemplate.hpp>

class FourNodeQuad: virtual public Element,
                            public QuadElementTemplate<double> 
{
 protected:
        int nn[4];
 public:
        FourNodeQuad(int*);

        Element *clone();

        void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);
        double getMassSensitivityWRTthickness(CoordSet&);
        double weight(CoordSet&, double *);
        double weightDerivativeWRTthickness(CoordSet&, double *, int=1); 

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *gs);
        void getGravityForceSensitivityWRTthickness(CoordSet&,double *gravity, int senMethod, Vector& f, int gravflg,
                                                    GeomState *gs);
        void             getVonMises (Vector &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
                                      double ylayer=0.0, double zlayer=0.0, int avgnum=0);
        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double>*, CoordSet &cs,
                                                Vector &elDisp, int strInd, int surface, int senMethod, double *ndTemps,
                                                int avgnum, double ylayer, double zlayer);

        void             getAllStress(FullM &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp,
                                      int strInd, int surface=0,
                                      double *ndTemps=0);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

        int getTopNumber();


        void computeDisp(CoordSet&cs, State &state, const InterpPoint &,
                         double*res, GeomState *gs);
        void getFlLoad(CoordSet &, const InterpPoint &,  double *flF,
                       double *resF, GeomState *gs=0);
        void getThermalForce(CoordSet &cs, Vector &ndTemps, 
                             Vector &ThermalForce, int glflag, GeomState *gs=0);
        int dim() { return 2; }

        void spy();
	
	// DEC
	PrioInfo examine(int sub, MultiFront *);

        Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};
#endif

