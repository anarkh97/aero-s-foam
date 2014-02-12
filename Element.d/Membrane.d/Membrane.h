#ifndef _MEMBRANE_H_
#define _MEMBRANE_H_

#include <Element.d/Element.h>
#include <Element.d/Membrane.d/MembraneElementTemplate.cpp>

class Membrane : public Element
#ifdef USE_EIGEN3
                 , public MembraneElementTemplate<double> 
#endif
{

        int nn[3];
  public:
        Membrane(int*);

        Element *clone();

        void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet& cs);
        double weight(CoordSet&, double *, int);
        double weightDerivativeWRTthickness(CoordSet&, double *, int, int=1);

        void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                                 GeomState *gs);

        void             getVonMises
                                (Vector &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface=0,
                                 double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void             getAllStress
                                (FullM &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface=0,
                                 double *ndTemps=0);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        int getTopNumber();

	// Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);

        bool hasRot() { return true; }

        int  getMassType() { return 0; } // lumped only

        Corotator* getCorotator(CoordSet &cs, double *kel, int fitAlg, int);
};
#endif
