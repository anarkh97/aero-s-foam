#ifndef _DISCRETEMASS6DOF_H_
#define _DISCRETEMASS6DOF_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

#include <Eigen/Core>

// Discrete mass+inertia with offset, no gravity force

class DiscreteMass6Dof : public Element, public Corotator
{
    int nn[1];
    Eigen::Matrix3d *C0;
    Eigen::Vector3d f0; // m*g;

  public:
    DiscreteMass6Dof(int*);
    ~DiscreteMass6Dof();

    void setFrame(EFrame *elemframe);
    void renum(int*);
    void renum(EleRenumMap&) override;
    int numNodes() const override { return 1; }
    int *nodes(int *) const override;
    int numDofs() const { return 6; }
    int* dofs(DofSetArray&, int* = 0);
    void markDofs(DofSetArray&);
    bool hasRot() { return true; }
    int getTopNumber() override { return 506; }

    FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
    int getMassType() { return 0; } // lumped only
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
    double getMass(CoordSet&);
    void getGravityForce(CoordSet&, double *g, Vector &f, int gravflg, GeomState *gs=0);
    void computePressureForce(CoordSet&, Vector& elPressureForce, GeomState *, int, double);

    // nonlinear functions
    Corotator* getCorotator(CoordSet&, double*, int, int) { return this; }
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &c0,
                          FullSquareMatrix &elk, double *f, double dt, double t);
    void getInternalForce(GeomState *, GeomState &curState, CoordSet &,
                          FullSquareMatrix&, double *f, double, double);
    double getElementEnergy(GeomState&, CoordSet&);
    bool useDefaultInertialStiffAndForce() { return false; }
    void getInertialStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0,
                                  FullSquareMatrix &elK, double *f, double dt, double t,
                                  double beta, double gamma, double alphaf, double alpham);
    void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                             int &nlflag) { nlflag = 2; }
    void getNLVonMises(Vector&, Vector& weight, GeomState &, CoordSet &, int);
};

#endif
