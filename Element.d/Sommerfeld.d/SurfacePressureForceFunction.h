#ifndef _SURFACEPRESSUREFORCEFUNCTION_H_
#define _SURFACEPRESSUREFORCEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Eigen/Geometry>

template<typename Scalar, template <typename S> class ShapeFunctionTemplate, typename QuadratureRule>
class SurfacePressureForceFunction
 : public VectorValuedFunction<(ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates+1)*ShapeFunctionTemplate<Scalar>::NumberOfValues,
                               (ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates+1)*ShapeFunctionTemplate<Scalar>::NumberOfValues,
                               Scalar,
                               (ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates+1)*ShapeFunctionTemplate<Scalar>::NumberOfValues+9,
                               2,
                               double>
{
  public:
    enum {
      NumberOfNodes = ShapeFunctionTemplate<Scalar>::NumberOfValues,
      NumberOfSurfaceDimensions = ShapeFunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      NumberOfDimensions = NumberOfSurfaceDimensions+1,
      NumberOfGeneralizedCoordinates = NumberOfNodes*NumberOfDimensions,
      NumberOfValues = NumberOfNodes*NumberOfDimensions,
      NumberOfScalarConstants = NumberOfNodes*NumberOfDimensions+9,
      NumberOfIntegerConstants = 2
    };

  private:
    Eigen::Matrix<double,NumberOfNodes,NumberOfDimensions> X;
    double p; // pressure
    int deg;  // quadrature rule degree
    BlastLoading::BlastData *conwep;

  public:
    SurfacePressureForceFunction(const Eigen::Array<double,NumberOfScalarConstants,1>& sconst,
                                 const Eigen::Array<int,NumberOfIntegerConstants,1>& iconst)
    {
      for(int inode = 0; inode < NumberOfNodes; ++inode) {
        X.row(inode) = Eigen::Map<Eigen::Matrix<double,NumberOfDimensions,1> >(const_cast<double*>(sconst.data())+NumberOfDimensions*inode);
      }
      deg = iconst[0];
      if(iconst[1] == 0) { // constant pressure
        p = sconst[NumberOfNodes*NumberOfDimensions];
        conwep = 0;
      }
      else {
        conwep = new BlastLoading::BlastData;
        p = sconst[NumberOfNodes*NumberOfDimensions];
        conwep->ExplosivePosition[0]    = sconst[NumberOfNodes*NumberOfDimensions+1];
        conwep->ExplosivePosition[1]    = sconst[NumberOfNodes*NumberOfDimensions+2];
        conwep->ExplosivePosition[2]    = sconst[NumberOfNodes*NumberOfDimensions+3];
        conwep->ExplosiveDetonationTime = sconst[NumberOfNodes*NumberOfDimensions+4];
        conwep->ExplosiveWeight         = sconst[NumberOfNodes*NumberOfDimensions+5];
        conwep->ScaleLength             = sconst[NumberOfNodes*NumberOfDimensions+6];
        conwep->ScaleTime               = sconst[NumberOfNodes*NumberOfDimensions+7];
        conwep->ScaleMass               = sconst[NumberOfNodes*NumberOfDimensions+8];
        conwep->BlastType = (iconst[1] == 1) ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst;
        conwep->ExplosiveWeightCubeRoot = pow(conwep->ExplosiveWeight,1.0/3.0);
      }
    }

    ~SurfacePressureForceFunction() {
      if(conwep) delete conwep;
    }

    Eigen::Matrix<Scalar,NumberOfValues,1> operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) const
    {
      // inputs:
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = x translation of node 2
      // q[4] = y translation of node 2
      // q[5] = z translation of node 2
      // etc...

      // return value: force vector
      Eigen::Matrix<Scalar,NumberOfValues,1> fext;
      fext.setZero();

      using std::abs;

      // Set current configuration
      Eigen::Matrix<Scalar,NumberOfNodes,NumberOfDimensions> x;
      for(int inode = 0; inode < NumberOfNodes; ++inode) {
        x.row(inode) = X.row(inode).template cast<Scalar>() + q.template segment<3>(inode*3).transpose();
      }

      // shape function and derivatives
      Eigen::Array<double, ShapeFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
      Eigen::Array<int, ShapeFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
      VectorValuedFunctionSpatialView<double, ShapeFunctionTemplate> S(sconst, iconst, 0);
      Eigen::AutoDiffJacobian<VectorValuedFunctionSpatialView<double, ShapeFunctionTemplate> > dN(S);

      // quadrature rule
      QuadratureRule c(deg);

      // local variables to be computed at each integration point
      Eigen::Matrix<double,NumberOfSurfaceDimensions,1> xi;          // abscissa
      double weight;                                                 // weight
      Eigen::Matrix<double,NumberOfNodes,1> N;                       // shape function values
      Eigen::Matrix<double,NumberOfNodes,NumberOfSurfaceDimensions> dNdXi;  // derivative of shape functions w.r.t. xi
      Eigen::Matrix<Scalar,NumberOfSurfaceDimensions,NumberOfDimensions> j;
      Eigen::Matrix<Scalar,NumberOfDimensions,1> normal;
      Eigen::Matrix<double,NumberOfSurfaceDimensions,NumberOfDimensions> J;
      Eigen::Matrix<double,NumberOfDimensions,1> Xip, Normal; // reference coordinates and normal
      double pressure = p;

      // Loop over the integration points
      for(int ip = 0; ip < c.getN(); ++ip) {

        // get the integration point abscissa and weight
        c.getAbscissaAndWeight(ip, xi, weight);

        // compute shape functions and derivatives
        dN(xi, &N, &dNdXi);

        // compute the surface normal
        j = dNdXi.template cast<Scalar>().transpose()*x;
        normal = j.row(0).cross(j.row(1));

        if(conwep) {
          J = dNdXi.transpose()*X;
          Normal = J.row(0).cross(J.row(1));
          Xip = X.transpose()*N;
          pressure = p + BlastLoading::ComputeGaussPointPressure(Xip.data(), Normal.data(), t, *conwep);
        }

        for(int i = 0; i < NumberOfDimensions; ++i)
          for(int k = 0; k < NumberOfNodes; ++k)
            fext[k*NumberOfDimensions+i] += Scalar(N(k)*weight*pressure)*normal(i);
      }

      return -fext;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
