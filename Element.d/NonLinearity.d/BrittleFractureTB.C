#include <Utils.d/NodeSpaceArray.h>
#include <Utils.d/pstress.h>
#include <cmath>

template<typename BaseMaterial>
int
BrittleFractureTB<BaseMaterial>::getNumStates()
{
  const int i = BaseMaterial::getNumStates();
  return i+1;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::initStates(double *state)
{
  BaseMaterial::initStates(state);
  const int i = BaseMaterial::getNumStates();
  state[i] = 0;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  // Note: this function is called for post-processing.
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(state[i] < Kf + tol) {
    BaseMaterial::getStress(_stress, _strain, state, temp);
  }
  else {
    Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
    stress->setZero();
  }
}

template<typename BaseMaterial>
void 
BrittleFractureTB<BaseMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &en, Tensor &_enp,
                                           double *staten, double *statenp, double temp, double dt)
{
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(staten[i] <  Kf + tol) {
    // compute elastic response
    BaseMaterial::integrate(_stress, _tm, en, _enp, staten, statenp, temp, dt);

    // check failure criteria
    std::vector<double> pvec(3);
    double svec[6];
    // Tensor_d0s2_Ss12 has xx, xy, xz, yy, yz, zz while pstress takes in xx, yy, zz, xy, yz, xz
    svec[0] = (*stress)[0];
    svec[1] = (*stress)[3];
    svec[2] = (*stress)[5];
    svec[3] = (*stress)[1];
    svec[4] = (*stress)[4];
    svec[5] = (*stress)[2];
    pstress(svec, pvec.data());
    double pmax = std::max(std::max(pvec[0],pvec[1]), pvec[2]); 
    if(pmax > maxprs) {
      statenp[i] = staten[i] + dt*pow(pmax-maxprs, exponent);
      if(statenp[i] >= Kf + tol) {
        stress->setZero();
        tm->setZero();
      }
    }
    else
      statenp[i] = staten[i];
  }
  else {
    statenp[i] = staten[i];
    stress->setZero();
    tm->setZero();
  }
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::integrate(Tensor *_stress, Tensor &en, Tensor &_enp,
                                           double *staten, double *statenp, double temp, double dt)
{
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(staten[i] <  Kf + tol) {
    // compute elastic response
    BaseMaterial::integrate(_stress, en, _enp, staten, statenp, temp, dt);

    // check failure criteria
    std::vector<double> pvec(3);
    double svec[6];
    // Tensor_d0s2_Ss12 has xx, xy, xz, yy, yz, zz while pstress takes in xx, yy, zz, xy, yz, xz
    svec[0] = (*stress)[0];
    svec[1] = (*stress)[3];
    svec[2] = (*stress)[5];
    svec[3] = (*stress)[1];
    svec[4] = (*stress)[4];
    svec[5] = (*stress)[2];
    pstress(svec, pvec.data());
    double pmax = std::max(std::max(pvec[0],pvec[1]), pvec[2]); 
    if(pmax > maxprs) {
      statenp[i] = staten[i] + dt*pow(pmax-maxprs, exponent);
      if(statenp[i] >= Kf + tol) {
        stress->setZero();
      }
    }
    else
      statenp[i] = staten[i];
  }
  else {
    statenp[i] = staten[i];
    stress->setZero();
  }
}

template<typename BaseMaterial>
double
BrittleFractureTB<BaseMaterial>::getDamage(double *statenp)
{
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();
  return (statenp[i] >= Kf + tol) ? 1.0 : 0.0;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::print(std::ostream &out) const
{
  BaseMaterial::print(out);
  out << " " << maxprs << " " << exponent << " " << Kf;
}
