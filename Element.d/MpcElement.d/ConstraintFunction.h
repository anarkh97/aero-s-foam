#ifndef _CONSTRAINTFUNCTION_H_
#define _CONSTRAINTFUNCTION_H_

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/AutoDiff>
//#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

// scalar valued rheonomic constraint, takes (x,t) as input
template<int _NumberOfGeneralizedCoordinates,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class RheonomicConstraintFunction
{
  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };
    virtual Scalar operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) const = 0;
};

// wrapper "Functor" to support automatic and numerical differentiation of constraint functions w.r.t spatial variables
template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class SpatialView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = 1
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;

  public:
    SpatialView(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                      ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const 
                      Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1>&
                      _iconst, Scalar _t = 0) 
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const 
    {
      ConstraintFunctionTemplate<T> F(sconst, iconst);
      f(0) = F(q, static_cast<T>(t));
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* f) const 
    { 
      return (*this)(q,*f);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class SpatialJacobian
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;

  public:
    SpatialJacobian(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                    ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                    Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                    _iconst, Scalar _t = 0)
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const
    {
      SpatialView<T, ConstraintFunctionTemplate> F(sconst, iconst, static_cast<T>(t));
      Eigen::AutoDiffJacobian<SpatialView<T, ConstraintFunctionTemplate> > J(F);
      Eigen::Matrix<T,1,ValuesAtCompileTime> ft;
      Eigen::Matrix<T,1,1> val;
      J(q,&val,&ft);
      for(int i=0; i<ValuesAtCompileTime; ++i) f(i) = ft(i);
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* f) const
    {
      return (*this)(q,*f);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// wrapper "Functor" to support automatic and numerical differentiation of constraint functions w.r.t t
template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class TemporalView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = 1,
      ValuesAtCompileTime = 1
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates, 1> q;

  public:
    TemporalView(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const 
                 Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1>&
                 _iconst, const Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::
                 NumberOfGeneralizedCoordinates, 1>& _q) 
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const 
    {
      ConstraintFunctionTemplate<T> F(sconst, iconst);
      f(0) = F(q.template cast<T>(),t(0));
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* f) const 
    { 
      return (*this)(t,*f);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class TemporalJacobian
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = 1,
      ValuesAtCompileTime = 1
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates, 1> q;

  public:
    TemporalJacobian(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                    ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                    Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                    _iconst, const Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::
                    NumberOfGeneralizedCoordinates, 1>& _q)
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& v) const
    {
      TemporalView<T, ConstraintFunctionTemplate> f(sconst, iconst, q.template cast<T>());
      Eigen::AutoDiffJacobian<TemporalView<T, ConstraintFunctionTemplate> > J(f);
      Eigen::Matrix<T,1,ValuesAtCompileTime> dfdt;
      Eigen::Matrix<T,1,1> val;
      J(t,&val,&dfdt);
      for(int i=0; i<ValuesAtCompileTime; ++i) v(i) = dfdt(i);
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* v) const
    {
      return (*this)(t,*v);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class TemporalViewOfSpatialJacobian
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = 1,
      ValuesAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates, 1> q;

  public:
    TemporalViewOfSpatialJacobian(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                    ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                    Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                    _iconst, const Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::
                    NumberOfGeneralizedCoordinates, 1>& _q)
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& J) const
    {
      SpatialJacobian<T,ConstraintFunctionTemplate> dfdq(sconst,iconst,t[0]);
      dfdq.template operator()<T>(q.template cast<T>(),J);
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* J) const
    {
      return (*this)(t,*J);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate>
class SpatialJacobianVelocityProduct
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = 1
    };

  private:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;
    Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates, 1> v;

  public:
    SpatialJacobianVelocityProduct(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                    ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                    Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                    _iconst, Scalar _t, Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::
                    NumberOfGeneralizedCoordinates, 1> _v)
     : sconst(_sconst), iconst(_iconst), t(_t), v(_v) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const
    {
      SpatialView<T, ConstraintFunctionTemplate> F(sconst, iconst, static_cast<T>(t));
      Eigen::AutoDiffJacobian<SpatialView<T, ConstraintFunctionTemplate> > dfdq(F);
      Eigen::Matrix<T,1,InputsAtCompileTime> J;
      Eigen::Matrix<T,1,1> val;
      dfdq(q,&val,&J);
      f[0] = J.dot(v.template cast<T>());
      return 1;
    }
    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* f) const
    {
      return (*this)(q,*f);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#ifdef USE_SACADO
#include <Sacado.hpp>

template<typename Functor> class SacadoReverseJacobian : public Functor
{
  typedef Sacado::RadVec::ADvar<typename Functor::Scalar> ActiveScalar;
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

  typedef Eigen::Matrix<ActiveScalar, InputType::SizeAtCompileTime, 1> ActiveInput;
  typedef Eigen::Matrix<ActiveScalar, Functor::ValueType::SizeAtCompileTime, 1> ActiveValue;

  template<typename T>
  int operator() (const Eigen::Matrix<T,InputType::SizeAtCompileTime,1>& x,
                  Eigen::Matrix<T,Functor::ValueType::SizeAtCompileTime,InputType::SizeAtCompileTime>& jac) const
  {

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

    return 0;
  }
};
#endif
#endif
#endif