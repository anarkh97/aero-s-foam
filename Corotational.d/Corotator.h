#ifndef _COROTATOR_H_
#define _COROTATOR_H_
#include <Utils.d/MyComplex.h>

class GeomState;
class TemperatureState;
class CoordSet;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
typedef GenFullM<DComplex> FullMC;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;

class Corotator {
  public:
    // DEFINITE
    virtual void getStiffAndForce(GeomState &, CoordSet &, 
                                  FullSquareMatrix &, double *, double, double)=0; 

    virtual void getDExternalForceDu(GeomState &geomState, CoordSet &cs,
                                     FullSquareMatrix &elK, double *locF){}

    virtual void getInternalForce(GeomState &geomState, CoordSet &cs,
                                  FullSquareMatrix &K, double *f, double dt, double t) 
    { getStiffAndForce(geomState, cs, K, f, dt, t); }

    virtual void getExternalForce(GeomState &, CoordSet &,  double *) {}

    // ONLY FOR BUCKLING (EIGEN)
    // ONLY NONLINEAR TERM OF K
    virtual void formGeometricStiffness(GeomState &, CoordSet &, 
                                        FullSquareMatrix &, double *); 

    virtual double* getOriginalStiffness();

    // ONLY FOR STRESSES
    virtual void extractDeformations(GeomState &geomState, CoordSet &cs, 
                                     double *vld, int &nlflag);
    virtual void extractDeformations(GeomState &geomState, CoordSet &cs,
                                     DComplex *vld, int &nlflag);

    virtual void getGlobalDisp(GeomState&, CoordSet&, Vector&) {}

    virtual void getNLVonMises(Vector&, Vector& weight,
                               GeomState &, CoordSet &, int);
    virtual void getNLVonMises(ComplexVector&, Vector& weight,
                               GeomState &, CoordSet &, int);

    virtual void getNLAllStress(FullM&, Vector&,
                                GeomState &, CoordSet &, int);
    virtual void getNLAllStress(FullMC&, Vector&,
                                GeomState &, CoordSet &, int);

    virtual double getElementEnergy(GeomState &, CoordSet &);
    virtual double getDissipatedEnergy(GeomState &, CoordSet &) { return 0.0; }
 
    // NOT USED
    virtual void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                        double *vlr);

    virtual void reBuildorigK(FullSquareMatrix&) {}

    // ONLY FOR MATERIAL NONLINEAR WITH INTERNAL VARIABLES
    virtual void getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &c0,
                                  FullSquareMatrix &elk, double *f, double dt, double t)
      { getStiffAndForce(curState, c0, elk, f, dt, t); }
    virtual void getInternalForce(GeomState *refState, GeomState &curState, CoordSet &c0,
                                  FullSquareMatrix &tmp, double *f, double dt, double t) 
      { getInternalForce(curState, c0, tmp, f, dt, t); }

    virtual void getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                               GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
                               double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
                               int avgnum = 0, int measure = -1)
      { getNLVonMises(stress, weight, curState, c0, strIndex); }

    virtual void getNLVonMises(ComplexVector& stress, Vector& weight, GeomState &curState,
                               GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
                               double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
                               int avgnum = 0, int measure = -1)
      { getNLVonMises(stress, weight, curState, c0, strIndex); }

    virtual void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
                                GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
                                double *ndTemps = 0, int measure = -1)
      { getNLAllStress(stress, weight, curState, c0, strInd); }

    virtual void getNLAllStress(FullMC &stress, Vector &weight, GeomState &curState,
                                GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
                                double *ndTemps = 0, int measure = -1)
      { getNLAllStress(stress, weight, curState, c0, strInd); }

    virtual void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0) {}

    virtual void getResidualCorrection(GeomState &gs, double *r) {}

    virtual void initMultipliers(GeomState& c1) {}
    virtual void updateMultipliers(GeomState& c1) {}
    virtual double getError() { return 0; }

    // ONLY USED FOR POSTPROCESSING TO OUTPUT STRESS OR STRAIN AT THE GAUSS POINTS
    // CURRENTLY ONLY SUPPORTED FOR MATERIAL NONLINEAR
    virtual int getNumGaussPoints() { return 0; }

    virtual ~Corotator() {/*TODO*/}
};

#endif
