#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Utils.d/dofset.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>

static HexahedralShapeFunction hexaSF;
static StructProp nullProp;

NLHexahedral::NLHexahedral(int *nd, bool isLinKin) {
  int i;
  for(i = 0; i < 8; ++i)
    n[i] = nd[i];
  linearKinematics = isLinKin;
  material = 0;
}

void
NLHexahedral::setProp(StructProp *p)
{
 material = new BilinPlasKinHardMat(p);
 prop = p;
}

void
NLHexahedral::setMaterial(Material *m)
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

void
NLHexahedral::getGaussPointAndWeight(int n, double *point, double &weight)
{
if (0==1){
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

  static double wi[3] =	 {
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

  static double wi[2] =	 {
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
NLHexahedral::getShapeFunction(){return &hexaSF;}

StrainEvaluator *
NLHexahedral::getStrainEvaluator(){
  if(linearKinematics)
    return &linearStrain;
  else
    return &greenLagrangeStrain;
}

Material *
NLHexahedral::getMaterial(){return material;}

int NLHexahedral::getTopNumber()
{
  if(linearKinematics)
    return 301;
  else
    return 302;
}
