#ifndef _SACADOREVERSEJACOBIAN_H_
#define _SACADOREVERSEJACOBIAN_H_

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>

#ifdef USE_SACADO
#include <Sacado.hpp>
#endif

/*
namespace Eigen {

namespace internal {
inline const Sacado::RadVec::ADvar<double>& conj(const Sacado::RadVec::ADvar<double>& x)  { return x; }
};

template<> struct NumTraits<Sacado::RadVec::ADvar<double> >
    : NumTraits<double>
{
  typedef Sacado::RadVec::ADvar<double> Real;
  typedef Sacado::RadVec::ADvar<double> NonInteger;
  typedef Sacado::RadVec::ADvar<double> Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};

template<int N> struct NumTraits<Sacado::Rad::ADvar<Sacado::Fad::SFad<double,N> > >
    : NumTraits<double>
{
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<double,N> > Real;
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<double,N> > NonInteger;
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<double,N> > Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};

}
*/

template<typename Functor, bool special_op=true> class SacadoReverseJacobian : public Functor
{
#ifdef USE_SACADO
  typedef Sacado::RadVec::ADvar<typename Functor::Scalar> ActiveScalar;
#endif
public:
  SacadoReverseJacobian() : Functor() {}
  SacadoReverseJacobian(const Functor& f) : Functor(f) {}

  // forward constructors
  template<typename T0>
  SacadoReverseJacobian(const T0& a0) : Functor(a0) {}
  template<typename T0, typename T1>
  SacadoReverseJacobian(const T0& a0, const T1& a1) : Functor(a0, a1) {}
  template<typename T0, typename T1, typename T2>
  SacadoReverseJacobian(const T0& a0, const T1& a1, const T2& a2) : Functor(a0, a1, a2) {}

  typedef typename Functor::Scalar Scalar;
  typedef typename Functor::InputType InputType;
  typedef typename Functor::JacobianType ValueType;
#ifdef USE_SACADO
  typedef Eigen::Matrix<ActiveScalar, InputType::SizeAtCompileTime, 1> ActiveInput;
  typedef Eigen::Matrix<ActiveScalar, Functor::ValueType::SizeAtCompileTime, 1> ActiveValue;
#endif

  template<typename T>
  int operator() (const Eigen::Matrix<T,InputType::SizeAtCompileTime,1>& x,
                  Eigen::Matrix<T,Functor::ValueType::SizeAtCompileTime,InputType::SizeAtCompileTime>& jac) const
  {
#ifdef USE_SACADO
// note: the 2nd condition is a hack to get code to build with buggy version of icpc 12
#if defined(_OPENMP) && (!defined(__INTEL_COMPILER) || __INTEL_COMPILER < 1200 || __INTEL_COMPILER > 1210)
    #pragma omp critical
    {
#endif
    ActiveInput ax = x.template cast<ActiveScalar>();
    ActiveValue av;

    Functor::operator()(ax, av);

    for (int i=0; i<jac.rows(); i++)
    {
      Sacado::RadVec::ADvar<typename Functor::Scalar>::Outvar_Gradcomp(av[i]);
      for (int j=0; j<jac.cols(); j++)
        jac.coeffRef(i,j) = ax[j].adj();
    }

    Sacado::RadVec::ADvar<typename Functor::Scalar>::aval_reset();
#if defined(_OPENMP) && (!defined(__INTEL_COMPILER) || __INTEL_COMPILER < 1200 || __INTEL_COMPILER > 1210)
    }
#endif
#else
    std::cerr << "warning: USE_SACADO is not defined\n";
#endif
    return 0;
  }
};

#ifdef USE_SACADO
template<typename Functor> class SacadoHessian : public Functor
{
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<typename Functor::Scalar,Functor::InputType::SizeAtCompileTime> > ActiveScalar;
public:

  SacadoHessian() : Functor() {}
  SacadoHessian(const Functor& f) : Functor(f) {}

  // forward constructors
  template<typename T0>
  SacadoHessian(const T0& a0) : Functor(a0) {}
  template<typename T0, typename T1>
  SacadoHessian(const T0& a0, const T1& a1) : Functor(a0, a1) {}
  template<typename T0, typename T1, typename T2>
  SacadoHessian(const T0& a0, const T1& a1, const T2& a2) : Functor(a0, a1, a2) {}

  typedef typename Functor::Scalar Scalar;
  typedef typename Functor::InputType InputType;
  typedef typename Functor::JacobianType ValueType;
  typedef typename Functor::HessianType JacobianType;

  typedef Eigen::Matrix<ActiveScalar, InputType::SizeAtCompileTime, 1> ActiveInput;
  //typedef Eigen::Matrix<ActiveScalar, Functor::ValueType::SizeAtCompileTime, 1> ActiveValue;
  typedef Eigen::Matrix<ActiveScalar, 1, Functor::ValueType::SizeAtCompileTime> ActiveValue;

  //int operator() (const InputType& x, ValueType& jac) const
  template<typename T>
  int operator() (const Eigen::Matrix<T,InputType::SizeAtCompileTime,1>& x, Eigen::Matrix<T,InputType::SizeAtCompileTime,InputType::SizeAtCompileTime>& hes) const
  {

    ActiveInput ax = x.template cast<ActiveScalar>();
    ActiveValue av;

    for (int j=0; j<hes.cols(); j++)
      ax[j] = Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime>(hes.cols(), j, x[j]);

    Functor::operator()(ax, av);

    Sacado::Rad::ADvar< Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime> >::Gradcomp();
    for (int i=0; i<hes.rows(); i++)
    {
      for (int j=0; j<hes.cols(); j++)
        hes.coeffRef(i,j) = ax[i].adj().dx(j);
    }

    Sacado::Rad::ADvar< Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime> >::aval_reset();

    return 0;
  }

};
#endif

#endif
#endif
