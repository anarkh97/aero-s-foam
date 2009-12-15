#ifndef PITA_PITANONLINDYNAM_H_
#define PITA_PITANONLINDYNAM_H_

#include <Problems.d/NonLinDynam.h>
#include <OOPita.d/DynamState.h>

/*#include <Pita.d/NiceTimer.h>*/

namespace Pita {

class PitaNonLinDynamic : public NonLinDynamic {
public:

  // TODO Timers
  // NiceTimerHandler pitaTimers;
  
  typedef Vector VecType;
  typedef double ScalarType;
  
  // Constructor
  PitaNonLinDynamic(Domain *);
   
  // Get initial displacement and velocity
  using NonLinDynamic::getInitState; // to avoid hiding
  int getInitState(DynamState &);
  int getInitSeedCount() const;
  int getInitSeed(DynamState &, int sliceRank);
 
  // Added Accessors
  int getKiter() const { return kiter; }
  int getJratio() const { return Jratio; }
  int getNumTSonCPU() const { return numTSonCPU; }
  int getNumTS() const { return numTS; }
  double getCoarseDt() const { return coarseDt; }
  double getCoarseDelta() const { return coarseDelta; }  
  const SparseMatrix * getStiffMatrix() const { return K; }
  int getBaseImprovementMethod() const { return baseImprovementMethod; }
 
  // Added methods
  void reBuildCoarse(GeomState & geomState, int iter = 0);
  void reBuildFine(GeomState & geomState, int iter = 0);
  void reBuildKonly();
  void zeroRotDofs(VecType &) const;
  double formRHSCoarseCorrector(Vector & inc_displac, Vector & velocity, Vector & acceleration, Vector & residual, Vector & rhs);
  void formRHSCoarsePredictor(Vector & velocity, Vector & residual, Vector & rhs, GeomState & geomState, double mid = 0.0);
  double energyNorm(const Vector & disp, const Vector & velo);
  double energyDot(const Vector & disp1, const Vector & velo1, const Vector & disp2, const Vector & velo2);

  // Output
  void openResidualFile();
  void pitaDynamOutput(int timeSliceRank, GeomState * geomState, Vector & velocity,
                       Vector & vp, double time, int step, Vector & force, Vector & aeroF);
  void openOutputFiles(int sliceRank);
  void closeOutputFiles(); 
  //void printNLPitaTimerFile(int CPUid);

  /*class PitaPostProcessor : public NLDynamPostProcessor
  {
  public:
    explicit PitaPostProcessor(PitaNonLinDynamic & probDesc);
    virtual ~PitaPostProcessor();
    int sliceRank() const { return sliceRank_; }
    void sliceRank(int rank);
    virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                                Vector& velocity, Vector& bkVelocity,
                                Vector& vp, Vector& bkVp, int step, int parity,
                                int aeroAlg) {};
    virtual void dynamOutput(GeomState * geomState, Vector & velocity, Vector & vp, double time, int step, Vector & force, Vector & aeroF) const;
  private:
    PitaNonLinDynamic & probDesc_;
    int sliceRank_;  
  };*/

  //virtual const NLDynamPostProcessor & defaultPostProcessor() const { return defaultPostProcessor_; }
  //virtual PitaPostProcessor & defaultPostProcessor() { return defaultPostProcessor_; }

protected:
  SparseMatrix *K;               // PITA requires to explicitely build the stiffness matrix
  int kiter, Jratio, numTSonCPU; // PITA main parameters from input file
  int numTS;                     // Total number of time-slices 
  double coarseDt, coarseDelta;  // Coarse time parameters
  int baseImprovementMethod;     // 0 = all seeds (global), 1 = increments only (local)

private:
  // Overloaded method, to build the stiffness matrix during NonLinDynamic::preProcess() 
  void buildOps(AllOps<double> &, double, double, double, Rbm *);
  //PitaPostProcessor defaultPostProcessor_;
};

inline void
PitaNonLinDynamic::reBuildFine(GeomState & geomState, int iter)
{
  reBuild(geomState, iter, getDelta());
}

inline void
PitaNonLinDynamic::reBuildCoarse(GeomState & geomState, int iter)
{
  reBuild(geomState, iter, getCoarseDelta());
}

inline double
PitaNonLinDynamic::formRHSCoarseCorrector(Vector & inc_displac, Vector & velocity, Vector & acceleration, Vector & residual, Vector & rhs)
{
  return formRHScorrector(inc_displac, velocity, acceleration, residual, rhs, getCoarseDelta());
}

inline void
PitaNonLinDynamic::formRHSCoarsePredictor(Vector & velocity, Vector & residual, Vector & rhs, GeomState & geomState, double mid)
{
  formRHSpredictor(velocity, residual, rhs, geomState, mid, getCoarseDelta());
}

/*inline void
PitaNonLinDynamic::PitaPostProcessor::dynamOutput(GeomState * geomState, Vector & velocity, Vector & vp, double time, int step, Vector & force, Vector & aeroF) const
{
  probDesc_.pitaDynamOutput(sliceRank_, geomState, velocity, vp, time, step, force, aeroF);
}*/

} // end namespace Pita

#endif
