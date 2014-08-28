#ifndef _MODAL_DESCR_H_
#define _MODAL_DESCR_H_

#include <Problems.d/ModalBase.h>
#include <Math.d/Vector.h>

#include <Hetero.d/FlExchange.h>

template<class VecType> class SysState;
class FlExchanger;
template <typename Scalar> struct AllSensitivities;

template <class Scalar>
class ModalDescr : public ModalBase {
/* problem descriptor for linear dynamics problems using the eigenmdoes of
     non-zero frequency as the solution space
*/
private:

  ModalOps modalOps;
  AllSensitivities<Scalar> *allSens;

public:

  ModalDescr() {}
  ModalDescr(Domain *d);

  void projectForce(Vector &fullV, Vector& modalV);
  void expand(const Vector &modalV, Vector& fullV);

  void preProcess();
  void preProcessSA() {  filePrint(stderr," ... ModalDescr::preProcessSA is not implemented\n");  exit(-1);  }
  void postProcessSA(ModalOps *,Vector &sol) {  filePrint(stderr," ... ModalDescr::postProcessSA is not implemented\n");  exit(-1);  }
  void processLastOutput();

  ModalDescr* getPostProcessor() { return this; }
  void getTimes(double &dt, double &tmax);
  void printFullNorm(Vector &){};
  void getInitState(SysState<Vector> &state);

  double getInitialForceNorm();

  int solVecInfo() { return numModes; }
  int masterSolVecInfo() { return numModes; }
  int bcInfo() {return 0;}

  void getConstForce(Vector &constF);
  void addConstForceSensitivity(Vector &gravityForceSen);
  void getSteadyStateParam(int &flag, int &min, int &max, double &tol);
  int getTimeIntegration() { return domain->solInfo().timeIntegration; }
  void getNewMarkParameters(double &beta, double &gamma,
    double &alphaf, double  &alpham);
  Domain *getDomain() { return domain; };
  AllSensitivities<Scalar> *getAllSensitivities() { return allSens; }
  void getInitialTime(int &tIndex, double &time);
  void getRayleighCoef(double &alpha) { alpha =  domain->solInfo().alphaDamp; }

  ModalOps* buildOps(double kcoef, double ccoef, double mcoef);

  void computeStabilityTimeStep(double &dt, ModalOps &){ /* leave blank */ }
  void reSolve(ModalOps*, Vector &, int, Vector &){ /* leave blank */ }
  void updateState(double, Vector &, Vector &){ /* leave blank */ }
  void push_forward(Vector &){ /* leave blank */ }
  void pull_back(Vector &){ /* leave blank */ }
  void getQuasiStaticParameters(double &maxVel, double &delta);
  int getFilterFlag() { return domain->solInfo().filterFlags; }
  void project(Vector &v) { /* leave blank */ }

  void getContactForce(Vector&, Vector&, Vector &extF, double, double, double) { extF.zero(); };
  void computeExtForce2(SysState<Vector>& state, Vector &extF,
                        Vector &constF, int tIndex, double time, Vector *aeroF = 0,
                        double gamma = 0.5, double alphaf = 0.5);
  void getAeroelasticForceSensitivity(int t_index, double t, Vector * aero_f=0, double gamma=0.5, double alphaf=0.5) {
      filePrint(stderr," ... ModalDescr::getAeroelasticForceSensitivity\n");  exit(-1);
  }

  void getInternalForce(Vector &d, Vector &f, double t, int tIndex);

  void printTimers(ModalOps *, double) { /* leave blank */}

  int getModeDecompFlag(){ return domain->solInfo().modeDecompFlag; }
  void modeDecomp(double time, int tIndex, Vector& dsp){ /* leave blank */ }

  void dynamOutput(int tIndex, double time, ModalOps &mOps, Vector &extF,
                   Vector *aeroF, SysState<Vector> &state);

  FlExchanger *flExchanger;

  void computeTimeInfo() {};
  int aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  int aeroSensitivityPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  int sendDisplacements(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  void sendNumParam(int numParam){ flExchanger->sendNumParam(numParam); }
  void sendRelativeResidual(double relres){ flExchanger->sendRelativeResidual(relres); }
  int cmdCom(int cmdFlag){ return flExchanger->cmdCom(cmdFlag); }
 
  void a5TimeLoopCheck(int& parity, double& t, double dt);
  void a5StatusRevise(int parity, SysState<Vector>& curState, SysState<Vector>& bkState);
  int getAeroAlg() { return domain->solInfo().aeroFlag; }
  void aeroSend(double time, Vector& d, Vector& v, Vector& a, Vector& v_p) { domain->aeroSend(d, v, a, v_p, bcx, vcx); }

  void thermoePreProcess(Vector& d_n, Vector& v_n, Vector& v_p) { domain->thermoePreProcess(); }
  int getThermoeFlag() { return domain->solInfo().thermoeFlag; }
  void thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p) { domain->thermohPreProcess(d_n, v_n, v_p, bcx); }
  int getThermohFlag() { return domain->solInfo().thermohFlag; }

  void aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p) { domain->aeroHeatPreProcess(d_n, v_n, v_p, bcx); }
  int getAeroheatFlag() { return domain->solInfo().aeroheatFlag; }
};

#ifdef _TEMPLATE_FIX_
  #include <Problems.d/ModalDescrTem.C>
#endif

#endif
