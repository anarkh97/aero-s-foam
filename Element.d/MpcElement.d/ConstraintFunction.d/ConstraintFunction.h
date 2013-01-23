#ifndef _CONSTRAINTFUNCTION_H_
#define _CONSTRAINTFUNCTION_H_

#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>
#include <Eigen/Sparse>

#include <Element.d/Function.d/SacadoReverseJacobian.h>

// scalar valued rheonomic constraint, takes (x,t) as input
template<int _NumberOfGeneralizedCoordinates,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class RheonomicConstraintFunction {

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

template<typename _Scalar, template <typename S> class ConstraintFunctionTemplate, bool special_op=true>
class ConstraintJacobian
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = ConstraintFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates
    };

  protected:
    Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                 ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;

  public:
    ConstraintJacobian(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
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
class PartialTimeDerivative
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
    PartialTimeDerivative(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
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
class TemporalViewOfConstraintJacobian
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
    TemporalViewOfConstraintJacobian(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
                    ConstraintFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                    Eigen::Array<int, ConstraintFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                    _iconst, const Eigen::Matrix<Scalar, ConstraintFunctionTemplate<Scalar>::
                    NumberOfGeneralizedCoordinates, 1>& _q)
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& J) const
    {
      ConstraintJacobian<T,ConstraintFunctionTemplate> dfdq(sconst,iconst,t[0]);
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
class ConstraintJacobianVelocityProduct
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
    ConstraintJacobianVelocityProduct(const Eigen::Array<typename ConstraintFunctionTemplate<Scalar>::ScalarConstantType,
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

#endif
#endif
