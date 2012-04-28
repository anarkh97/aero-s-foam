#ifndef _SHAPE_FUNCTIONS_H_
#define _SHAPE_FUNCTIONS_H_

enum RegionType { LineSegment, Triangle, Quadrilateral, Tetrahedron, Wedge, Pyramid, Hexahedron };
enum InterpolationType { LagrangePolynomial, HermiteBirkhoffType1, HermiteBirkhoffType2, HermiteBirkhoffType3,
                         HermiteBirkhoffType1LS, HermiteBirkhoffType2LS, HermiteBirkhoffType3LS,
                         HermiteBirkhoffType1WLS, HermiteBirkhoffType2WLS, HermiteBirkhoffType3WLS };

template<RegionType Region>
struct region_traits { };

template<>
struct region_traits<LineSegment> {
  enum {
    NumberOfDimensions = 1,
    NumberOfVertices = 2
  };
};

template<>
struct region_traits<Triangle> {
  enum {
    NumberOfDimensions = 2,
    NumberOfVertices = 3
  };
};

template<>
struct region_traits<Quadrilateral> {
  enum {
    NumberOfDimensions = 2,
    NumberOfVertices = 4
  };
};

template<>
struct region_traits<Tetrahedron> {
  enum {
    NumberOfDimensions = 3,
    NumberOfVertices = 4
  };
};

template<>
struct region_traits<Wedge> {
  enum {
    NumberOfDimensions = 3,
    NumberOfVertices = 6 
  };
};

template<>
struct region_traits<Hexahedron> {
  enum {
    NumberOfDimensions = 3,
    NumberOfVertices = 8 
  };
};


#ifdef USE_EIGEN3
#include <Eigen/Core>

template<RegionType Region, int NY, typename _Scalar = double, InterpolationType Interpolation = LagrangePolynomial>
class ShapeFunctions
{
protected:
  int m_inputs, m_values;

public:
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = region_traits<Region>::NumberOfDimensions,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
  typedef int HessianType;

  ShapeFunctions() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  ShapeFunctions(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  template<typename T>
  int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& x, Eigen::Matrix<T,ValuesAtCompileTime,1> *y) const;
};
#endif

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Meta.d/ShapeFunctions.C>
#endif
#endif
