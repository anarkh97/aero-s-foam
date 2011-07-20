#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Utils.d/dofset.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>

void
HexahedralShapeFunction::getLocalDerivatives(Tensor *_localDerivatives, double xi[3])
{
  Tensor_d1s2_sparse *localDerivatives  = static_cast<Tensor_d1s2_sparse *>(_localDerivatives);

  static double localNodeCoordinates[24];
  localNodeCoordinates[0] = -1;
  localNodeCoordinates[1] = -1;
  localNodeCoordinates[2] = -1;
  localNodeCoordinates[3] = 1;
  localNodeCoordinates[4] = -1;
  localNodeCoordinates[5] = -1;
  localNodeCoordinates[6] = 1;
  localNodeCoordinates[7] = 1;
  localNodeCoordinates[8] = -1;
  localNodeCoordinates[9] = -1;
  localNodeCoordinates[10] = 1;
  localNodeCoordinates[11] = -1;
  localNodeCoordinates[12] = -1;
  localNodeCoordinates[13] = -1;
  localNodeCoordinates[14] = 1;
  localNodeCoordinates[15] = 1;
  localNodeCoordinates[16] = -1;
  localNodeCoordinates[17] = 1;
  localNodeCoordinates[18] = 1;
  localNodeCoordinates[19] = 1;
  localNodeCoordinates[20] = 1;
  localNodeCoordinates[21] = -1;
  localNodeCoordinates[22] = 1;
  localNodeCoordinates[23] = 1;

  // localDerivatives(k,i,j) = at node k, the derivative of shape function i w.r.t xi[j]
  for(int k=0; k<8; k++)
    for(int i=0; i<3; i++) {
      (*localDerivatives)[k][3*i  ]=(1.0/8.0)*(localNodeCoordinates[3*k  ]*(1+localNodeCoordinates[3*k+1]*xi[1])*(1+localNodeCoordinates[3*k+2]*xi[2]));
      (*localDerivatives)[k][3*i+1]=(1.0/8.0)*(localNodeCoordinates[3*k+1]*(1+localNodeCoordinates[3*k  ]*xi[0])*(1+localNodeCoordinates[3*k+2]*xi[2]));
      (*localDerivatives)[k][3*i+2]=(1.0/8.0)*(localNodeCoordinates[3*k+2]*(1+localNodeCoordinates[3*k  ]*xi[0])*(1+localNodeCoordinates[3*k+1]*xi[1]));
    }
}


static HexahedralShapeFunction hexaSF;
static StructProp nullProp;

NLHexahedral::NLHexahedral(int *nd, int _strainMeasure)
{
  int i;
  for(i = 0; i < 8; ++i)
    n[i] = nd[i];
  strainMeasure = _strainMeasure;
  material = 0;
}

void
NLHexahedral::setProp(StructProp *p)
{
  cerr << "here in NLHexahedral::setProp\n";
  material = new BilinPlasKinHardMat(p);
  prop = p;
}

void
NLHexahedral::setMaterial(NLMaterial *m)
{
  material = m;
  prop = &nullProp;
}

int
NLHexahedral::getNumGaussPoints()
{
  if (0==1)
    return 3*3*3;
  else
    return 2*2*2;
}

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;
extern DeformationGradient deformationGradient;

void
NLHexahedral::getGaussPointAndWeight(int n, double *point, double &weight)
{
  if (0==1) {
    int i, j, k;
    i = n%3;
    n /= 3;
    j = n%3;
    n /= 3;
    k = n;

    static double xi[3] = {
  	-0.7745966692,
	0.,
	0.7745966692
    };

    static double wi[3] = {
	5.0/9.0,
	8.0/9.0,
	5.0/9.0
    };
    weight = wi[i]*wi[j]*wi[k];
    point[0] = xi[i];
    point[1] = xi[j];
    point[2] = xi[k];
  }
  else {

    int i, j, k;
    i = n%2;
    n /= 2;
    j = n%2;
    n /= 2;
    k = n;

    static double xi[2] = {
  	-0.5773502692,
	 0.5773502692
    };

    static double wi[2] = {
	1,
	1
    };
    weight = wi[i]*wi[j]*wi[k];
    point[0] = xi[i];
    point[1] = xi[j];
    point[2] = xi[k];
  }
}

void
NLHexahedral::renum(int *table)
{
  int i;
  for(i =0; i < 8; ++i)
    n[i] = table[n[i]];
}

void
NLHexahedral::markDofs(DofSetArray &dsa)
{
  dsa.mark(n, 8, DofSet::XYZdisp);
}

int *
NLHexahedral::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[24];
  int i;
  for(i = 0; i < 8; ++i)
    dsa.number(n[i], DofSet::XYZdisp, p+3*i);
  return p;
}

int*
NLHexahedral::nodes(int *p)
{
  if(p == 0) p = new int[8];
  int i;
  for(i = 0; i < 8; ++i)
    p[i] = n[i];
  return p;
}

ShapeFunction *
NLHexahedral::getShapeFunction()
{
  return &hexaSF;
}

StrainEvaluator *
NLHexahedral::getStrainEvaluator()
{
  switch(strainMeasure) {
    case 0: return &linearStrain;
    case 1: return &greenLagrangeStrain;
    case 2: return &deformationGradient;
  }
  return material->getStrainEvaluator(); // default
}

NLMaterial *
NLHexahedral::getMaterial()
{
  return material;
}

int NLHexahedral::getTopNumber()
{
  return 117;
}

