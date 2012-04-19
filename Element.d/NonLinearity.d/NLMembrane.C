#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Utils.d/dofset.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Math.d/TTensor.h>
#include <Corotational.d/utilities.h>
#include <Corotational.d/GeomState.h>

template <int n>
void
LinearStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor  &e,
                        typename TwoDTensorTypes<n>::GradUTensor &gradU)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
}


template <int n>
void
LinearStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                              typename TwoDTensorTypes<n>::BTensor &B, 
			      typename TwoDTensorTypes<n>::DBTensor &DB,
			      typename TwoDTensorTypes<n>::GradUTensor &gradU, 
			      typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0];
    B[k][1] = dgradUdqk[k][1][1];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]);
  }
}

template <int n>
void
GLStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor &e,
                    typename TwoDTensorTypes<n>::GradUTensor &gradU)
{
  // 1/2*((I_32+gradU)^T(I+gradU)-I_32)
  e[0] = gradU[0][0] +
         0.5*(gradU[0][0]*gradU[0][0]+gradU[0][1]*gradU[0][1]+gradU[0][2]*gradU[0][2]);
  e[1] = gradU[1][1] +
         0.5*(gradU[1][0]*gradU[1][0]+gradU[1][1]*gradU[1][1]+gradU[1][2]*gradU[1][2]);
  e[2] = 0.5*(gradU[0][1]+gradU[1][0] +
         gradU[1][0]*gradU[0][0]+gradU[1][1]*gradU[0][1]+gradU[1][2]*gradU[0][2]);
}

template <int n>
void
GLStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                          typename TwoDTensorTypes<n>::BTensor &B, 
			  typename TwoDTensorTypes<n>::DBTensor &DB,
			  typename TwoDTensorTypes<n>::GradUTensor &gradU, 
			  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  // 1/2*(I+gradU)^T(I+gradU)-I
  e[0] = gradU[0][0] +
         0.5*(gradU[0][0]*gradU[0][0]+gradU[0][1]*gradU[0][1]+gradU[0][2]*gradU[0][2]);
  e[1] = gradU[1][1] +
         0.5*(gradU[1][0]*gradU[1][0]+gradU[1][1]*gradU[1][1]+gradU[1][2]*gradU[1][2]);
  e[2] = 0.5*(gradU[0][1]+gradU[1][0] +
           gradU[1][0]*gradU[0][0]+gradU[1][1]*gradU[0][1]+gradU[1][2]*gradU[0][2]);
  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0]
        + gradU[0][0]*dgradUdqk[k][0][0]
	+ gradU[0][1]*dgradUdqk[k][0][1]
	+ gradU[0][2]*dgradUdqk[k][0][2];
    B[k][1] = dgradUdqk[k][1][1]
        + gradU[1][0]*dgradUdqk[k][1][0]
	+ gradU[1][1]*dgradUdqk[k][1][1]
	+ gradU[1][2]*dgradUdqk[k][1][2];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]
        + gradU[0][0]*dgradUdqk[k][1][0]
	+ gradU[0][1]*dgradUdqk[k][1][1]
	+ gradU[0][2]*dgradUdqk[k][1][2]
	+ gradU[1][0]*dgradUdqk[k][0][0]
	+ gradU[1][1]*dgradUdqk[k][0][1]
	+ gradU[1][2]*dgradUdqk[k][0][2]);
  }
  int l;
  for(k = 0; k < n; ++k) {
    for(l=0; l < n; ++l) {
      DB[k][l][0] = dgradUdqk[k][0][0]*dgradUdqk[l][0][0]
                  + dgradUdqk[k][0][1]*dgradUdqk[l][0][1]
                  + dgradUdqk[k][0][2]*dgradUdqk[l][0][2];
      DB[k][l][1] = dgradUdqk[k][1][0]*dgradUdqk[l][1][0]
                  + dgradUdqk[k][1][1]*dgradUdqk[l][1][1]
                  + dgradUdqk[k][1][2]*dgradUdqk[l][1][2];
      DB[k][l][2] = 0.5*( dgradUdqk[k][0][0]*dgradUdqk[l][1][0]
                  + dgradUdqk[k][0][1]*dgradUdqk[l][1][1]
                  + dgradUdqk[k][0][2]*dgradUdqk[l][1][2]
		  + dgradUdqk[l][0][0]*dgradUdqk[k][1][0]
                  + dgradUdqk[l][0][1]*dgradUdqk[k][1][1]
                  + dgradUdqk[l][0][2]*dgradUdqk[k][1][2]);
    }
  }

}

//typedef SymTensor<double,2> Stress2D;
//typedef SimpleTensor<Stress2D, 9> Stress2DDeriv9;
void
TriMembraneShapeFunct::getGlobalGrads(Grad2D &gradU, Grad2DDeriv9 &dGradUdqk,
                                      double *_jac,
                                      Node *nodes, double xi[3], Vector &disp)
{
  int i,j,k;
  // First obtain the frame X, in which xi[i] runs along X[i].
  double d[2][3] = { { nodes[1].x-nodes[0].x, 
                       nodes[1].y-nodes[0].y,
                       nodes[1].z-nodes[0].z },
                     { nodes[2].x-nodes[0].x, 
                       nodes[2].y-nodes[0].y,
                       nodes[2].z-nodes[0].z },
		   };
  //fprintf(stderr, "Dirs %e %e %e   %e %e %e\n", d[0][0], d[0][1],d[0][2], d[1][0],d[1][1],
  //d[1][2]);
  double X[3][3];
  X[2][0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  X[2][1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  X[2][2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
 
  double l1inv = 1.0/sqrt(d[0][0]*d[0][0]+d[0][1]*d[0][1]+d[0][2]*d[0][2]);
  double l3inv = 1.0/sqrt(X[2][0]*X[2][0]+X[2][1]*X[2][1]+X[2][2]*X[2][2]);
  //fprintf(stderr, "linvs: %e %e\n", l1inv, l3inv);
  for(i = 0; i < 3; ++i) {
    X[0][i] = l1inv*d[0][i];
    X[2][i] *= l3inv;
  }
  X[1][0] = X[2][1]*X[0][2] - X[2][2]*X[0][1];
  X[1][1] = X[2][2]*X[0][0] - X[2][0]*X[0][2];
  X[1][2] = X[2][0]*X[0][1] - X[2][1]*X[0][0];
 
  // fprintf(stderr, "Axis %e %e %e   %e %e %e %e %e %e\n", X[0][0], X[0][1], X[0][2], 
  //   X[1][0], X[1][1], X[1][2], X[2][0], X[2][1], X[2][2]);
 
  // This is dX_i/dxi_j
  double dXdxi[2][2] = { { d[0][0]*X[0][0]+d[0][1]*X[0][1]+d[0][2]*X[0][2],
                           d[0][0]*X[1][0]+d[0][1]*X[1][1]+d[0][2]*X[1][2] },
                         { d[1][0]*X[0][0]+d[1][1]*X[0][1]+d[1][2]*X[0][2],
                           d[1][0]*X[1][0]+d[1][1]*X[1][1]+d[1][2]*X[1][2] }
		       };
  double &jac = *_jac;
  jac = dXdxi[0][0]*dXdxi[1][1]-dXdxi[0][1]*dXdxi[1][0];
  double jacInv = 1.0/jac;
  double dxidX[2][2] = { { jacInv*dXdxi[0][0], -jacInv*dXdxi[1][0] },
                         { -jacInv*dXdxi[0][1], jacInv*dXdxi[1][1] }
                       };
  // Create the local displacement vectors.
  double Utilde[3][3] = 
    { { X[0][0]*disp[0] + X[0][1]*disp[1] + X[0][2]*disp[2],
        X[1][0]*disp[0] + X[1][1]*disp[1] + X[1][2]*disp[2],
	X[2][0]*disp[0] + X[2][1]*disp[1] + X[2][2]*disp[2]
      },
      { X[0][0]*disp[3] + X[0][1]*disp[4] + X[0][2]*disp[5],
        X[1][0]*disp[3] + X[1][1]*disp[4] + X[1][2]*disp[5],
	X[2][0]*disp[3] + X[2][1]*disp[4] + X[2][2]*disp[5]
      },
      { X[0][0]*disp[6] + X[0][1]*disp[7] + X[0][2]*disp[8],
        X[1][0]*disp[6] + X[1][1]*disp[7] + X[1][2]*disp[8],
	X[2][0]*disp[6] + X[2][1]*disp[7] + X[2][2]*disp[8]
      } };
  // Derivatives of the scalar shape functions:
  double dNdX[3][2] = {  { -(dxidX[0][0] + dxidX[1][0]), 
                           -(dxidX[0][1] + dxidX[1][1]) },
                         { dxidX[0][0], dxidX[0][1] },
                         { dxidX[1][0], dxidX[1][1] }
		      };
  // fprintf(stderr, "dNdX: %e %e  %e %e   %e %e\n", 
  //     dNdX[0][0], dNdX[0][1], dNdX[1][0], dNdX[1][1], dNdX[2][0], dNdX[2][1]);
  for(i = 0; i < 2; ++i) {
    gradU[i][0] = dNdX[0][i]*Utilde[0][0]
                + dNdX[1][i]*Utilde[1][0]
                + dNdX[2][i]*Utilde[2][0];
    gradU[i][1] = dNdX[0][i]*Utilde[0][1]
                + dNdX[1][i]*Utilde[1][1]
                + dNdX[2][i]*Utilde[2][1];
    gradU[i][2] = dNdX[0][i]*Utilde[0][2]
                + dNdX[1][i]*Utilde[1][2]
                + dNdX[2][i]*Utilde[2][2];
  }

  for(k = 0; k < 9; ++k) {
    int dir, nd;
    dir = k%3;
    nd = k/3;
    for(i = 0; i < 2; ++i)
      for(j = 0; j < 3; ++j) {
        dGradUdqk[k][i][j] = dNdX[nd][i]*X[j][dir];
      }
  }
  // The volume of the element is only half of the jacobian.
  jac *=0.5;
}

void
TriMembraneShapeFunct::getGradU(Grad2D &gradU,
                       Node *nodes, double xi[3], Vector &disp)
{
  int i;
  // First obtain the frame X, in which xi[i] runs along X[i].
  double d[2][3] = { { nodes[1].x-nodes[0].x, 
                       nodes[1].y-nodes[0].y,
                       nodes[1].z-nodes[0].z },
                     { nodes[2].x-nodes[0].x, 
                       nodes[2].y-nodes[0].y,
                       nodes[2].z-nodes[0].z },
		   };
  // fprintf(stderr, "Dirs %e %e %e   %e %e %e\n", d[0][0], d[0][1],d[0][2], d[1][0],d[1][1],
  // d[1][2]);
  double X[3][3];
  X[2][0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  X[2][1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  X[2][2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
 
  double l1inv = 1.0/sqrt(d[0][0]*d[0][0]+d[0][1]*d[0][1]+d[0][2]*d[0][2]);
  double l3inv = 1.0/sqrt(X[2][0]*X[2][0]+X[2][1]*X[2][1]+X[2][2]*X[2][2]);
  // fprintf(stderr, "linvs: %e %e\n", l1inv, l3inv);
  for(i = 0; i < 3; ++i) {
    X[0][i] = l1inv*d[0][i];
    X[2][i] *= l3inv;
  }
  X[1][0] = X[2][1]*X[0][2] - X[2][2]*X[0][1];
  X[1][1] = X[2][2]*X[0][0] - X[2][0]*X[0][2];
  X[1][2] = X[2][0]*X[0][1] - X[2][1]*X[0][0];
 
  // fprintf(stderr, "Axis %e %e %e   %e %e %e %e %e %e\n", X[0][0], X[0][1], X[0][2], 
  //     X[1][0], X[1][1], X[1][2], X[2][0], X[2][1], X[2][2]);
 
  // This is dX_i/dxi_j
  double dXdxi[2][2] = { { d[0][0]*X[0][0]+d[0][1]*X[0][1]+d[0][2]*X[0][2],
                           d[0][0]*X[1][0]+d[0][1]*X[1][1]+d[0][2]*X[1][2] },
                         { d[1][0]*X[0][0]+d[1][1]*X[0][1]+d[1][2]*X[0][2],
                           d[1][0]*X[1][0]+d[1][1]*X[1][1]+d[1][2]*X[1][2] }
		       };
  double jac;
  jac = dXdxi[0][0]*dXdxi[1][1]-dXdxi[0][1]*dXdxi[1][0];
  double jacInv = 1.0/jac;
  double dxidX[2][2] = { { jacInv*dXdxi[0][0], -jacInv*dXdxi[1][0] },
                         { -jacInv*dXdxi[0][1], jacInv*dXdxi[1][1] }
                       };
  // Create the local displacement vectors.
  double Utilde[3][3] = 
    { { X[0][0]*disp[0] + X[0][1]*disp[1] + X[0][2]*disp[2],
        X[1][0]*disp[0] + X[1][1]*disp[1] + X[1][2]*disp[2],
	X[2][0]*disp[0] + X[2][1]*disp[1] + X[2][2]*disp[2]
      },
      { X[0][0]*disp[3] + X[0][1]*disp[4] + X[0][2]*disp[5],
        X[1][0]*disp[3] + X[1][1]*disp[4] + X[1][2]*disp[5],
	X[2][0]*disp[3] + X[2][1]*disp[4] + X[2][2]*disp[5]
      },
      { X[0][0]*disp[6] + X[0][1]*disp[7] + X[0][2]*disp[8],
        X[1][0]*disp[6] + X[1][1]*disp[7] + X[1][2]*disp[8],
	X[2][0]*disp[6] + X[2][1]*disp[7] + X[2][2]*disp[8]
      } };
  // Derivatives of the scalar shape functions:
  double dNdX[3][2] = {  { -(dxidX[0][0] + dxidX[1][0]), 
                           -(dxidX[0][1] + dxidX[1][1]) },
                         { dxidX[0][0], dxidX[0][1] },
                         { dxidX[1][0], dxidX[1][1] }
		      };
  // fprintf(stderr, "dNdX: %e %e  %e %e   %e %e\n", 
  //     dNdX[0][0], dNdX[0][1], dNdX[1][0], dNdX[1][1], dNdX[2][0], dNdX[2][1]);
  for(i = 0; i < 2; ++i) {
    gradU[i][0] = dNdX[0][i]*Utilde[0][0]
                + dNdX[1][i]*Utilde[1][0]
                + dNdX[2][i]*Utilde[2][0];
    gradU[i][1] = dNdX[0][i]*Utilde[0][1]
                + dNdX[1][i]*Utilde[1][1]
                + dNdX[2][i]*Utilde[2][1];
    gradU[i][2] = dNdX[0][i]*Utilde[0][2]
                + dNdX[1][i]*Utilde[1][2]
                + dNdX[2][i]*Utilde[2][2];
  }
}

static TriMembraneShapeFunct shpFct;

NLMembrane::NLMembrane(int *nd)
 : material(NULL)
{
  for(int i = 0; i < 3; ++i)
    n[i] = nd[i];
}

NLMembrane::~NLMembrane()
{
  if(material && useDefaultMaterial) delete material;
}

int
NLMembrane::getNumGaussPoints()
{
  return 1;
  //return 3;
}

void
NLMembrane::getGaussPointAndWeight(int n, double *point, double &weight)
{
  const double third = 1.0/3.0;
  point[0] = third;
  point[1] = third;
  point[2] = 0.0;
  weight = 1.0;
/*
  double w_save[3] = {
    0.33333333333333333333,
    0.33333333333333333333,
    0.33333333333333333333 };
  double xy_save[2*3] = {
    0.66666666666666666667,  0.16666666666666666667,
    0.16666666666666666667,  0.66666666666666666667,
    0.16666666666666666667,  0.16666666666666666667 };

  point[0] = xy_save[2*n+0];
  point[1] = xy_save[2*n+1];
  point[2] = 0.0;
  weight = w_save[n];
*/
}

void
NLMembrane::renum(int *table)
{
  n[0] = table[n[0]];
  n[1] = table[n[1]];
  n[2] = table[n[2]];
}

int*
NLMembrane::nodes(int *nd)
{
  if(nd == 0) nd = new int[3];
  nd[0] = n[0];
  nd[1] = n[1];
  nd[2] = n[2];
  return nd;
}

int *
NLMembrane::dofs(DofSetArray &dsa, int*df)
{
  if(df == 0) df = new int[9];
  dsa.number(n[0], DofSet::XYZdisp, df);
  dsa.number(n[1], DofSet::XYZdisp, df+3);
  dsa.number(n[2], DofSet::XYZdisp, df+6);

  return df;
}

LinearStrain2D<9> linStrain2D;
GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
NLMembrane::getGenStrainEvaluator()
{
  return material->getGenStrainEvaluator();
}

NLMaterial *
NLMembrane::getMaterial()
{
  return material;
}

GenShapeFunction< TwoDTensorTypes<9> > *
NLMembrane::getShapeFunction()
{
 return &shpFct;
}

void
NLMembrane::markDofs(DofSetArray &dsa)
{
  dsa.mark(n, 3, DofSet::XYZdisp);
}

void
NLMembrane::setMaterial(NLMaterial *m)
{
  material = m;
}

int
NLMembrane::numInternalNodes()
{
  // this function is called after setMaterial
  useDefaultMaterial = (material == NULL);
  if(useDefaultMaterial) material = new ElaLinIsoMat2D(prop);
  return 0;
}

void
NLMembrane::computePressureForce(CoordSet& cs, Vector& force,
                                 GeomState *geomState, int cflg, double)
{
  Node nodes[3];
  double gs[9];
  for(int i = 0; i < 3; ++i) { 
    nodes[i] = *cs[n[i]];
    gs[3*i  ] = (*geomState)[n[i]].x;
    gs[3*i+1] = (*geomState)[n[i]].y;
    gs[3*i+2] = (*geomState)[n[i]].z;
  }

  double d[2][3] = { { nodes[1].x+gs[3]-nodes[0].x-gs[0], 
                       nodes[1].y+gs[4]-nodes[0].y-gs[1],
                       nodes[1].z+gs[5]-nodes[0].z-gs[2] },
                     { nodes[2].x+gs[6]-nodes[0].x-gs[0], 
                       nodes[2].y+gs[7]-nodes[0].y-gs[1],
                       nodes[2].z+gs[8]-nodes[0].z-gs[2] },
		   };
  double n[3];
  n[0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  n[1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  n[2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
  double p = getPressure();
  for(int i=0; i < 3; ++i)
    force[i] = force[i+3]=force[i+6] = 1.0/6.0*p*n[i];
}

#include <Corotational.d/PhantomCorotator.h>
#include <Corotational.d/MatNLCorotator.h>
Corotator*
NLMembrane::getCorotator(CoordSet &, double *, int , int)
{
  if(prop == NULL) return new PhantomCorotator();
  else return new MatNLCorotator(this, false);
}

FullSquareMatrix
NLMembrane::stiffness(CoordSet& cs, double *k, int flg)
{
  if(prop == NULL) { 
    FullSquareMatrix ret(9,k);
    ret.zero();
    return ret;
  }
  else {
    return GenGaussIntgElement<TwoDTensorTypes<9> >::stiffness(cs,k,flg);
  }
}

#include <Eigen/Core>

FullSquareMatrix
NLMembrane::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix ret(9,mel);
  if(prop == NULL) { ret.zero(); return ret; }

  if(cmflg) { // consistent mass matrix

    double mass = getMass(cs);
    Eigen::Map<Eigen::Matrix<double,9,9> > M(mel);
    M << 2, 0, 0, 1, 0, 0, 1, 0, 0,
         0, 2, 0, 0, 1, 0, 0, 1, 0,
         0, 0, 2, 0, 0, 1, 0, 0, 1,
         1, 0, 0, 2, 0, 0, 1, 0, 0,
         0, 1, 0, 0, 2, 0, 0, 1, 0,
         0, 0, 1, 0, 0, 2, 0, 0, 1,
         1, 0, 0, 1, 0, 0, 2, 0, 0,
         0, 1, 0, 0, 1, 0, 0, 2, 0,
         0, 0, 1, 0, 0, 1, 0, 0, 2;
     M *= mass/12;
  }
  else { // lumped mass matrix

    double mass = getMass(cs);
    double massPerNode = mass/3.0;

    ret.zero();
    for(int i = 0; i < 9; ++i)
      ret[i][i] = massPerNode;
  }

  return ret;
}

double
NLMembrane::getMass(CoordSet& cs)
{
  if(prop == NULL) return 0;

  Node &nd1 = cs.getNode(n[0]);
  Node &nd2 = cs.getNode(n[1]);
  Node &nd3 = cs.getNode(n[2]);

  double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

  v1[0] = r3[0] - r1[0];
  v1[1] = r3[1] - r1[1];
  v1[2] = r3[2] - r1[2];

  v2[0] = r2[0] - r1[0];
  v2[1] = r2[1] - r1[1];
  v2[2] = r2[2] - r1[2];

  crossprod(v1, v2, v3);

  double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
  double density = prop->rho;
  double t       = prop->eh;

  double mass = area*t*density;

  return mass;
}

void
NLMembrane::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                            Vector& gravityForce, int gravflg, GeomState *geomState)
{
  // TODO
  gravityForce.zero();
}

void
NLMembrane::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                        Vector& elDisp, int strInd, int, double *ndTemps,
                        double ylayer, double zlayer, int avgnum)
{
  // TODO
  stress.zero();
  weight.zero();
}

void
NLMembrane::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                         Vector& elDisp, int strInd, int, double *ndTemps)
{
  // TODO
  stress.zero();
  weight.zero();
}

#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
void
NLMembrane::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip,
                        double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(n[0], xyz[0], xyz[0]+3);
  state.getDV(n[1], xyz[1], xyz[1]+3);
  state.getDV(n[2], xyz[2], xyz[2]+3);

  for(int j = 0; j < 6; ++j)
    res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
}

void
NLMembrane::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                      double *resF, GeomState *gs)
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]   = (1.0-gp[0]-gp[1]) * flF[i];
    resF[3+i] = gp[0] * flF[i];
    resF[6+i] = gp[1] * flF[i];
  }
}

// Four node membrane comprising two three node membranes, 3 dof per node
NLMembrane4::NLMembrane4(int *nodenums)
{
  int i,j,k;
  nn = new int[4];
  for(i=0; i<4; ++i) nn[i] = nodenums[i];

  nSubElems = 2;
  subElems = new Element * [2];
  subElemNodes = new int * [2];
  subElemDofs = new int * [2];

  subElemNodes[0] = new int[3];
  subElemNodes[0][0] = 0; subElemNodes[0][1] = 1; subElemNodes[0][2] = 3;

  subElemNodes[1] = new int[3];
  subElemNodes[1][0] = 2; subElemNodes[1][1] = 3; subElemNodes[1][2] = 1;

  for(i=0; i<2; ++i) {
    int tmp[3];
    subElemDofs[i] = new int[9];
    for(j=0; j<3; ++j) {
      int nij = subElemNodes[i][j];
      tmp[j] = nodenums[nij]; // global node numbers
      for(k=0;k<3;++k) {
        subElemDofs[i][3*j+k] = 3*nij+k;
      }
    }
    subElems[i] = new NLMembrane(tmp);
  }
  nnodes = 4;
  ndofs = 12;
}

int
NLMembrane4::getTopNumber()
{
  return 102;
}

void
NLMembrane4::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  subElems[i]->computeDisp(cs, state, subip, res);
}

void
NLMembrane4::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                       double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  double subres[9];
  subElems[i]->getFlLoad(cs, subip, flF, subres);

  int j;
  for(j=0; j<12; ++j) res[j] = 0.0;
  for(j=0; j<9; ++j) res[subElemDofs[i][j]] = subres[j];
}

