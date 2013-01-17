#ifndef _VECTORVALUEDFUNCTION_H_
#define _VECTORVALUEDFUNCTION_H_

#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>
#include <Eigen/Sparse>

#include <Element.d/Function.d/SacadoReverseJacobian.h>

// Vector-valued temporospatial function, takes (q,t) as input where q is a Vector (spatial) and t is a scalar (temporal)
template<int _NumberOfGeneralizedCoordinates,
         int _NumberOfValues,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class VectorValuedFunction {

  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfValues                 = _NumberOfValues,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };

    virtual Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) const = 0;
};

// Vector-valued temporospatial function, takes (q,t) as input where q is a Matrix (spatial) and t is a scalar (temporal)
template<int InputNumberOfRows, int InputNumberOfColumns,
         int _NumberOfValues,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class VectorValuedFunctionOfAMatrix 
: public VectorValuedFunction<InputNumberOfRows*InputNumberOfColumns,
                              _NumberOfValues, 
                              _Scalar,
                              _NumberOfScalarConstants,
                              _NumberOfIntegerConstants,
                              _ScalarConstantType>
{
  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { NumberOfGeneralizedCoordinates = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues                 = _NumberOfValues,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };

    virtual Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) const = 0;

    Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) const
    {
      Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns> qmat;
      for(int i=0; i<InputNumberOfRows; ++i)
        for(int j=0; j<InputNumberOfColumns; ++j)
          qmat(i,j) = q[i*InputNumberOfColumns + j];

      return this->operator() (qmat, t);
    }
};

// Vector-valued temporospatial function, takes (q,t) as input where q is a scalar (spatial) and t is a scalar (temporal)
template<int _NumberOfValues,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class VectorValuedFunctionOfAScalar
: public VectorValuedFunction<1,
                              _NumberOfValues,                          
                              _Scalar,
                              _NumberOfScalarConstants,
                              _NumberOfIntegerConstants,
                              _ScalarConstantType>
{
  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { NumberOfGeneralizedCoordinates = 1,
           NumberOfValues                 = _NumberOfValues,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };

    virtual Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Scalar q, Scalar t) const = 0;

    Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) const
    {
      return this->operator() (q[0], t);
    }
};


// wrapper "Functor" to support automatic and numerical differentiation of vector valued temporospatial functions w.r.t spatial inputs, q
template<typename _Scalar, template <typename S> class VectorValuedFunctionTemplate>
class VectorValuedFunctionSpatialView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = VectorValuedFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = VectorValuedFunctionTemplate<Scalar>::NumberOfValues
    };

  private:
    Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                 VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;

  public:
    VectorValuedFunctionSpatialView(const Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const 
                Eigen::Array<int, VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1>&
                _iconst, Scalar _t = 0) 
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const 
    {
      VectorValuedFunctionTemplate<T> F(sconst, iconst);
      f = F(q, static_cast<T>(t));
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
    typedef Eigen::Array<Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime>,InputsAtCompileTime,1> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class VectorValuedFunctionTemplate, bool special_op=true>
class VectorValuedFunctionJacobian
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = VectorValuedFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = VectorValuedFunctionTemplate<Scalar>::NumberOfValues*
                            VectorValuedFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates
    };

  protected:
    Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                 VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    Eigen::Array<int,
                 VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> iconst;
    Scalar t;

  public:
    VectorValuedFunctionJacobian(const Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                       VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                       Eigen::Array<int, VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> &
                       _iconst, Scalar _t = 0)
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& f) const
    {
      VectorValuedFunctionSpatialView<T, VectorValuedFunctionTemplate> F(sconst, iconst, static_cast<T>(t));
      Eigen::AutoDiffJacobian<VectorValuedFunctionSpatialView<T, VectorValuedFunctionTemplate> > dFdq(F);
      Eigen::Matrix<T,VectorValuedFunctionTemplate<Scalar>::NumberOfValues,
                      VectorValuedFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> J;
      Eigen::Matrix<T,VectorValuedFunctionTemplate<Scalar>::NumberOfValues,1> val;
      dFdq(q,&val,&J);
      for(int i=0; i<ValuesAtCompileTime; ++i) f(i) = J.data()[i];
/*
      // sacado reverse autodiff: this is much slower. also it doesn't work with cross product function in Eigen
      VectorValuedFunctionSpatialView<T, VectorValuedFunctionTemplate> F(sconst, iconst, static_cast<T>(t));
      SacadoReverseJacobian<VectorValuedFunctionSpatialView<T, VectorValuedFunctionTemplate> > dFdq(F);
      Eigen::Matrix<T,VectorValuedFunctionTemplate<Scalar>::NumberOfValues,
                      VectorValuedFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> J;
      Eigen::Matrix<T,VectorValuedFunctionTemplate<Scalar>::NumberOfValues,1> val;
      dFdq(q,J);
      for(int i=0; i<ValuesAtCompileTime; ++i) f(i) = J.data()[i];
*/
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
