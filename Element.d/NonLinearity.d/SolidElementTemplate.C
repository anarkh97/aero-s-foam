#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/dofset.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#endif
#include <Element.d/Meta.d/ShapeFunctions.h>

template<RegionType Region, int NumberOfNodes>
void
ShapeFunctionTemplate<Region,NumberOfNodes>::getLocalDerivatives(Tensor *_localDerivatives, double _xi[3])
{
#ifdef USE_EIGEN3
  Eigen::AutoDiffJacobian<ShapeFunctions<Region,NumberOfNodes,double> > dN;
  Eigen::Matrix<double,3,1> xi;
  Eigen::Matrix<double,NumberOfNodes,1> N;      // shape function values
  Eigen::Matrix<double,NumberOfNodes,3> dNdXi;  // derivative of shape functions w.r.t. xi

  xi << _xi[0], _xi[1], _xi[2];

  dN(xi, &N, &dNdXi);

  Tensor_d1s2_sparse *localDerivatives  = static_cast<Tensor_d1s2_sparse *>(_localDerivatives);
  for(int k = 0; k < NumberOfNodes; ++k)
    for(int i = 0; i < 3; ++i) {
      (*localDerivatives)[k][3*i  ] = dNdXi(k,0);
      (*localDerivatives)[k][3*i+1] = dNdXi(k,1);
      (*localDerivatives)[k][3*i+2] = dNdXi(k,2);
    }
#else
  std::cerr << "ERROR: Implementation of ShapeFunctionTemplate::getLocalDerivatives in file" << std::endl
            << "       Element.d/NonLinearity.d/SolidElementTemplate.C required Eigen3 library.\n";
  exit(-1);
#endif
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::SolidElementTemplate(int *nd) 
 : material(0)
{
  for(int i = 0; i < NumberOfNodes; ++i)
    n[i] = nd[i];
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::setMaterial(NLMaterial *m)
{
  material = m;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
int
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getNumGaussPoints()
{
  return NumIntgPts;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
const ShapeFunctionTemplate<Region,NumberOfNodes> SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::shapeFunction;

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
ShapeFunction *
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getShapeFunction()
{
  return const_cast<ShapeFunctionTemplate<Region,NumberOfNodes>* >(&shapeFunction);
}

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;
extern DeformationGradient deformationGradient;

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
StrainEvaluator *
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getStrainEvaluator()
{
  return material->getStrainEvaluator();
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
NLMaterial *
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getMaterial()
{
  return material;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getNodeRefCoords(double (*_nodeRefCoords)[3])
{
  for(int i = 0; i < NumberOfNodes; ++i)
    for(int j = 0; j < 3; ++j)
      _nodeRefCoords[i][j] = nodeRefCoords[i][j];
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::getLocalNodalCoords(int i, double *coords)
{
  for(int j = 0; j < 3; ++j)
    coords[j] = nodeRefCoords[i][j];
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
int
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::numNodes()
{
  return NumberOfNodes;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
int 
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::numDofs()
{
  return NumberOfNodes*3;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::renum(int *table)
{
  for(int i = 0; i < NumberOfNodes; ++i)
    n[i] = table[n[i]];
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
void
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::markDofs(DofSetArray &dsa)
{
  dsa.mark(n, NumberOfNodes, DofSet::XYZdisp);
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
int *
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[3*NumberOfNodes];
  for(int i = 0; i < NumberOfNodes; ++i)
    dsa.number(n[i], DofSet::XYZdisp, p+3*i);
  return p;
}

template<RegionType Region, int NumberOfNodes, int NumIntgPts>
int*
SolidElementTemplate<Region,NumberOfNodes,NumIntgPts>::nodes(int *p)
{
  if(p == 0) p = new int[NumberOfNodes];
  for(int i = 0; i < NumberOfNodes; ++i)
    p[i] = n[i];
  return p;
}

/*
typedef SolidElementTemplate<Hexahedron,8, 8 > NLHexahedral8;
typedef SolidElementTemplate<Hexahedron,20,27> NLHexahedral20;
typedef SolidElementTemplate<Hexahedron,32,64> NLHexahedral32;

typedef SolidElementTemplate<Wedge,6, 6 > NLPentahedral6;
typedef SolidElementTemplate<Wedge,15,9 > NLPentahedral15;
typedef SolidElementTemplate<Wedge,26,18> NLPentahedral26;

typedef SolidElementTemplate<Tetrahedron,4 ,1 > NLTetrahedral4;
typedef SolidElementTemplate<Tetrahedron,10,15> NLTetrahedral10;
*/
