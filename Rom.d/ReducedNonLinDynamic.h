#ifndef ROM_REDUCEDNONLINDYNAMIC_H
#define ROM_REDUCEDNONLINDYNAMIC_H

#include "VecBasis.h"

#include <Math.d/Vector.h>

#include <memory>
#include <vector>

class Domain;
class GeomState;
class Corotator;

template <typename Scalar> class GenFullSquareMatrix;
template <typename Scalar> class GenSparseMatrix;

namespace Rom {

class ReducedNonLinDynamic {
public:
  explicit ReducedNonLinDynamic(Domain *d);
  ~ReducedNonLinDynamic();

  void preProcess();
  void computeTimeInfo();
  
  void getInitialTime(int &initTimeIndex, double &initTime) const;
  int getInitState(GenVector<double> &d, GenVector<double> &v, GenVector<double> &a, GenVector<double> &v_p) const;
  
  void readRestartFile(GenVector<double> &d, GenVector<double> &v, GenVector<double> &a, GenVector<double> &v_p, GeomState &);

  int solVecInfo() const;
  int sysVecInfo() const;
  int elemVecInfo() const;

  double getDt() const;
  double getDelta() const;
  int getMaxStep() const;

  void getNewmarkParameters(double &beta, double &gamma, double &alpha_f, double &alpha_m) const;
  
  double getTolerance() const;
  int getMaxit() const;
  
  int getAeroAlg() const;
  int getThermoeFlag() const;
  int getThermohFlag() const;
  int getAeroheatFlag() const;

  class Solver; 
  Solver *getSolver();
  const Solver *getSolver() const;
  
  GeomState *createGeomState() const;
  GeomState *copyGeomState(const GeomState *geomState) const;
  void updateStates(GeomState *refState, GeomState &geomState);
  
  void updatePrescribedDisplacement(GeomState *geomState) const;
  
  void getConstForce(GenVector<double> &constantForce) const;
  void getExternalForce(GenVector<double> &externalForce, const GenVector<double> &constantForce, int tIndex, double time,
                        const GeomState *geomState, GenVector<double> &elementInternalForce, const GenVector<double> &aeroF, double delta) const;
  
  // Form element stiffness matrices and return the residual force = external - internal forces
  double getStiffAndForce(GeomState &geomState, GenVector<double> &residual, GenVector<double> &elementInternalForce, double midtime, GeomState *refState);

  // Assemble new dynamic stiffness matrix
  void reBuild(GeomState& geomState, int iter, double delta, double t);

  void formRHSinitializer(GenVector<double> &fext, GenVector<double> &velocity, GenVector<double> &elementInternalForce, GeomState &geomState, GenVector<double> &rhs); // dummy
  void formRHSinitializer(GenVector<double> &fext, GenVector<double> &velocity, GenVector<double> &elementInternalForce, GeomState &geomState, GenVector<double> &rhs, GeomState *refState);
  double formRHScorrector(GenVector<double> &inc_displacement, GenVector<double> &velocity, GenVector<double> &acceleration, GenVector<double> &residual, GenVector<double> &rhs);

  int checkConvergence(int iter, double rN, GenVector<double> &residual, GenVector<double> &dv, double time);

  void dynamOutput(GeomState *geomState, GenVector<double> &velocity, GenVector<double> &vp,
                   double time, int timestep, GenVector<double> &force, GenVector<double> &aeroF, GenVector<double> &acceleration,
                   GeomState *refState);
  void processLastOutput();
  
  void dynamCommToFluid(GeomState *geomState, GeomState *bkGeomState,
                        GenVector<double> &velocity, GenVector<double> &bkVelocity,
                        GenVector<double> &vp, GenVector<double> &bkVp,
                        int step, int parity, int aeroAlg);
  
  void printTimers(double timeLoop);

  double getResidualNorm(GenVector<double> &res) { return res.norm(); }

private:
  int fullSolVecInfo() const;
  bool hasDamping() const;
  
  // Add to residual the internal and follower forces
  void getStiffAndReducedForceFromDomain(Vector &residual,
                                         GeomState &geomState, double time, GeomState *refState,
                                         Vector &elementInternalForce);

  double getIncrementTolerance() const;

  void zeroRotDofs(GenVector<double> &) const;

  Domain *domain_;

  int timeStepCount_;

  GenVecBasis<double, GenVector> reducedBasis_;

  std::vector<Corotator *> corotators_;
  GenFullSquareMatrix<double> *melArray_;
  GenFullSquareMatrix<double> *kelArray_;
  
  std::auto_ptr<GenSparseMatrix<double> > fullMatrix_;

  std::vector<double> vcx_;

  std::auto_ptr<GenFullSquareMatrix<double> > reducedMass_;
  std::auto_ptr<GenFullSquareMatrix<double> > reducedDamping_;
  std::auto_ptr<GenFullSquareMatrix<double> > reducedJacobian_;

  double firstResidualNorm_, firstIncrementNorm_, firstEnergyNorm_, secondResidualNorm_;

  std::auto_ptr<GenFullSquareMatrix<double> > zeroRotMetric_;

public:
  // Auxiliary classes
  class Solver {
  public:
    void reSolve(GenVector<double> &);

  private:
    const GenFullSquareMatrix<double> &matrix_;

    explicit Solver(const GenFullSquareMatrix<double> &);

    friend class ReducedNonLinDynamic;

    // Disallow copy & assignment
    Solver(const Solver &);
    Solver& operator=(const Solver &);
  };

  class Updater;

private:
  Solver reducedSolver_;

  friend class Updater;

  // Disallow copy & assigment
  ReducedNonLinDynamic(const ReducedNonLinDynamic &);
  ReducedNonLinDynamic &operator=(const ReducedNonLinDynamic &);
};


class ReducedNonLinDynamic::Updater {
public:
  typedef ReducedNonLinDynamic ProbDescr;
  typedef GenVector<double> VecType;
  typedef GeomState GeomType;

  typedef VecType StateIncr;
  typedef GeomType RefState;
  
  static RefState *initRef(const GeomType *u);

  static StateIncr *initInc(const GeomType *, const VecType *v);

  static void copyState(const GeomType *gn, RefState *gp);

  static void zeroInc(StateIncr *du);

  static void get_inc_displacement(ProbDescr *pbd,
                                   GeomType *geomState, StateIncr &du, GeomType *refState,
                                   bool zeroRot); 

  static double integrate(ProbDescr *pbd, RefState *refState, GeomType *geomState,
	                     	  StateIncr *du, VecType &residual,
                          VecType &elementInternalForce, VecType &gRes, VecType& vel_n,
                          VecType &accel, double midTime);

  static void midpointIntegrate(ProbDescr *pbd,
                                VecType &velN, double delta, GeomType *refState,
                                GeomType *geomState,
                                StateIncr *, VecType &,
                                VecType &, VecType &, VecType &acceleration,
                                bool zeroRot);

  static void updateIncr(StateIncr *du, VecType &ddu);

  static double formRHScorrector(ProbDescr *pbd, VecType &inc_displac, VecType &vel_n,
                                 VecType &accel, VecType &residual, VecType &rhs, GeomType *geomState, double delta);

  static void copyTo(RefState *, GeomType *, GeomType *, StateIncr *, VecType &, VecType &, VecType &, VecType &,
                     RefState *, GeomType *, GeomType *, StateIncr *, VecType &, VecType &, VecType &, VecType &);
};

} // end namespace Rom

#endif /* ROM_REDUCEDNONLINDYNAMIC_H */
