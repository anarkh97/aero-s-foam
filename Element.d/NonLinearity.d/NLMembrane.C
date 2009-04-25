#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Utils.d/dofset.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Math.d/TTensor.h>

template <int n>
void LinearStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor  &e,
 typename TwoDTensorTypes<n>::GradUTensor &gradU)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
}


template <int n>
void LinearStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
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
void GLStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor  &e,
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
void GLStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
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


class TriMembraneShapeFunct : public GenShapeFunction< TwoDTensorTypes<9> >{
  public:
    TriMembraneShapeFunct() : GenShapeFunction< TwoDTensorTypes<9> >() {}
    void getGlobalGrads(Grad2D &gradU, Grad2DDeriv9 &dGradUdqk, double *_jac,
                        Node *nodes, double xi[3], Vector &disp);
    void getGradU(Grad2D &gradU,
                        Node *nodes, double xi[3], Vector &disp);	
};

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
{ int i;
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

NLMembrane::NLMembrane(int *nd, bool isLinKin) {
  int i;
  for(i = 0; i < 3; ++i)
    n[i] = nd[i];
  linearKinematics = isLinKin;
  material = 0;
}


void
NLMembrane::setProp(StructProp *p)
{
 //material = new BilinPlasKinHardMat(p);
 prop = p;
}

int
NLMembrane::getNumGaussPoints()
{
 return 1;
}

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;

void
NLMembrane::getGaussPointAndWeight(int n, double *point, double &weight)
{
 const double third = 1.0/3.0;
 point[0] = third;
 point[1] = third;
 point[2] = 0.0;
 weight = 1.0;
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

StrainEvaluator *
NLMembrane::getStrainEvaluator()
{
  return 0;
}

static LinearStrain2D<9> linStrain;
static GLStrain2D<9> glStrain;

GenStrainEvaluator<TwoDTensorTypes<9> > *
NLMembrane::getGenStrainEvaluator()
{
//  fprintf(stderr, linearKinematics ? "Linear Strain\n" : "Non Linear Strain\n");
  return linearKinematics ? (GenStrainEvaluator<TwoDTensorTypes<9> > *)&linStrain 
            : (GenStrainEvaluator<TwoDTensorTypes<9> > *)&glStrain;
}

Material *
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
        dsa.mark(n, 3,  DofSet::XYZdisp);
}

void
NLMembrane::setMaterial(Material *m)
{
 //fprintf(stderr, "Setting material to %p\n", m);
 material = m;
}

void
NLMembrane::computePressureForce(Node *nodes,Vector& force,
                                          double *gs, int cflg)
{
 
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
