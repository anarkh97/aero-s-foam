#ifndef _SPACEJACOBIAN_H_
#define _SPACEJACOBIAN_H_

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>
#include <iostream>

#include <Element.d/Function.d/SacadoReverseJacobian.h>

namespace Simo {

template<typename A, typename B>
struct assign_coherent_impl {
  static void run(const A& a, B& b) { b = a; }
};

template<typename A, typename B>
void assign_coherent(const A& a, B& b)
{
  assign_coherent_impl<A,B>::run(a, b);
}

template<typename Scalar, int Options, int MaxRows, int MaxCols>
struct assign_coherent_impl<Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols>, Scalar> {
  typedef Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> A;
  typedef Scalar B;
  static void run(const A& a, B& b) { b = a[0]; }
};

template<typename Scalar, int Options, int MaxRows, int MaxCols>
struct assign_coherent_impl<Scalar, Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> > {
  typedef Scalar A;
  typedef Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> B;
  static void run(const A& a, B& b) { b[0] = a; }
};

template<typename Scalar, int A_Rows, int A_Cols, int A_Options, int A_MaxRows, int A_MaxCols,
                          int B_Rows, int B_Cols, int B_Options, int B_MaxRows, int B_MaxCols>
struct assign_coherent_impl<Eigen::Matrix<Scalar, A_Rows, A_Cols, A_Options, A_MaxRows, A_MaxCols>, 
                            Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols> > {
  typedef Eigen::Matrix<Scalar, A_Rows, A_Cols, A_Options, A_MaxRows, A_MaxCols> A;
  typedef Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols> B;
  static void run(const A& a, B& b) { b = Eigen::Map<B>(const_cast<Scalar*>(a.data()),a.rows(),a.cols()); }
};

// wrapper "Functor" to support automatic and numerical differentiation of spatio-temporal
// matrix valued function of a matrix w.r.t spatial coordinates, q
template<typename _Scalar, template <typename S> class FunctionTemplate>
class SpatialView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = FunctionTemplate<Scalar>::NumberOfValues
    };

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;
    Scalar t;

  public:
    SpatialView(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const 
                Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>&
                _iconst, Scalar _t = 0) 
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& _q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& y) const 
    {
      Eigen::Matrix<T,FunctionTemplate<T>::InputNumberOfRows,FunctionTemplate<T>::InputNumberOfColumns> q;
      assign_coherent(_q, q);

      // evaluate y = f(q)
      FunctionTemplate<T> f(sconst, iconst);
      assign_coherent(f(q,static_cast<T>(t)), y);

      return 1;
    }

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* y) const 
    { 
      return (*this)(q,*y);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    typedef Eigen::Array<Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime>,InputsAtCompileTime,1> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// wrapper "Functor" to support automatic and numerical differentiation of spatio-temporal
// matrix valued function of a matrix w.r.t the temporal coordinate, t
template<typename _Scalar, template <typename S> class FunctionTemplate>
class TemporalView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = 1,
      ValuesAtCompileTime = FunctionTemplate<Scalar>::NumberOfValues
    };

  private:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;
    const Eigen::Matrix<Scalar, FunctionTemplate<Scalar>::InputNumberOfRows, FunctionTemplate<Scalar>::InputNumberOfColumns>& q;

  public:
    TemporalView(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType, FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst,
                 const Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst,
                 const Eigen::Matrix<Scalar, FunctionTemplate<Scalar>::InputNumberOfRows, FunctionTemplate<Scalar>::InputNumberOfColumns >& _q) 
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& y) const 
    {
      // evaluate y = f(t)
      FunctionTemplate<T> f(sconst, iconst);
      assign_coherent(f(q.template cast<T>(), t[0]), y);

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

template<typename _Scalar, template <typename S> class FunctionTemplate, int method=0>
class SpaceJacobian
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows              = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns           = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants        = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants       = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues                 = FunctionTemplate<Scalar>::NumberOfValues*NumberOfGeneralizedCoordinates
    };
    typedef Eigen::Array<typename FunctionTemplate<Scalar>::ReturnType,InputNumberOfColumns,InputNumberOfRows> ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    SpaceJacobian(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                     FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                     Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                     _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::SpaceJacobianType
      operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

        typename FunctionTemplate<Scalar>::SpaceJacobianType ret;

        SpatialView<Scalar, FunctionTemplate> f(sconst, iconst, t);
        Eigen::AutoDiffJacobian<SpatialView<Scalar, FunctionTemplate> > dfdq(f);
        Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,
                      FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> J;
        Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,1> y;
        Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,1> _q;
        assign_coherent(q, _q);
        dfdq(_q, &y, &J);
        assign_coherent(J, ret);

        return ret;
      }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
#endif
