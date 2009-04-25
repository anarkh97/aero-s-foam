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
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;

class Corotator {
  public:
    // DEFINITE
    virtual void getStiffAndForce(GeomState &, CoordSet &, 
                                  FullSquareMatrix &, double *) = 0; 

    // ONLY FOR BUCKLING (EIGEN)
    // ONLY NONLINEAR TERM OF K
    virtual void formGeometricStiffness(GeomState &, CoordSet &, 
                                        FullSquareMatrix &, double *) { }; 

    virtual double* getOriginalStiffness() { return 0; };

    // ONLY FOR STRESSES
    virtual void extractDeformations(GeomState &geomState, CoordSet &cs, 
                                     double *vld, int &nlflag) { };
    virtual void extractDeformations(GeomState &geomState, CoordSet &cs,
                                     DComplex *vld, int &nlflag) { };

    virtual void getNLVonMises(Vector&, Vector& weight,
                               GeomState &, CoordSet &, int) { };
    virtual void getNLVonMises(ComplexVector&, Vector& weight,
                               GeomState &, CoordSet &, int) { };

    virtual void getNLAllStress(FullM&, Vector&,
                                GeomState &, CoordSet &, int) { };

    virtual double getElementEnergy(GeomState &, CoordSet &) { return 0.0; }; 
 
    // NOT USED
    virtual void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                        double *vlr) { };

    virtual ~Corotator() {/*TODO*/}
};

#endif
