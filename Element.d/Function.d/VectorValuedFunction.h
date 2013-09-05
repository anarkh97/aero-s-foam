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
    enum { InputNumberOfRows              = _NumberOfGeneralizedCoordinates,
           InputNumberOfColumns           = 1,
           NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfValues                 = _NumberOfValues,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };
    typedef Eigen::Matrix<Scalar,NumberOfValues,1> ReturnType;
    typedef Eigen::Matrix<Scalar,NumberOfValues,NumberOfGeneralizedCoordinates> SpaceJacobianType;

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

#endif
#endif
