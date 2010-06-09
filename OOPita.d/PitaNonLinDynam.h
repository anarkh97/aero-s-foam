#ifndef PITA_PITANONLINDYNAM_H_
#define PITA_PITANONLINDYNAM_H_

#include <Problems.d/NonLinDynam.h>
#include <OOPita.d/DynamState.h>

/*#include <Pita.d/NiceTimer.h>*/

namespace Pita {

class PitaNonLinDynamic : public NonLinDynamic {
public:
  typedef Vector VecType;
  typedef double ScalarType;

  // Overriden
  virtual void preProcess();

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
  void reBuildKonly();
  void zeroRotDofs(VecType &) const;
  double energyNorm(const Vector & disp, const Vector & velo);
  double energyDot(const Vector & disp1, const Vector & velo1, const Vector & disp2, const Vector & velo2);

  // Output
  void openResidualFile();
  void pitaDynamOutput(int timeSliceRank, GeomState * geomState, Vector & velocity,
                       Vector & vp, double time, int step, Vector & force, Vector & aeroF);
  void openOutputFiles(int sliceRank);
  void closeOutputFiles(); 

protected:
  SparseMatrix *K;               // PITA requires to explicitely build the stiffness matrix
  int kiter, Jratio, numTSonCPU; // PITA main parameters from input file
  int numTS;                     // Total number of time-slices 
  double coarseDt, coarseDelta;  // Coarse time parameters
  int baseImprovementMethod;     // 0 = all seeds (global), 1 = increments only (local)

private:
  // Overloaded method, to build the stiffness matrix during NonLinDynamic::preProcess() 
  void buildOps(AllOps<double> &, double, double, double, Rbm *);
};

} // end namespace Pita

#endif
