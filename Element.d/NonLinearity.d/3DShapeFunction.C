#include <stdio.h>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>

Tensor *
ShapeFunction::getGradUInstance()
{
  return (new Tensor_d0s2());
}

Tensor *
ShapeFunction::getDgradUDqkInstance()
{
  return (new Tensor_d1s2_sparse(numdofs));
}

void
ShapeFunction::getGlobalGrads(Tensor *_gradU, Tensor *_dgradUdqk, double *jac,
                              Node *nodes, double xi[3], Vector &disp)
{
  Tensor_d1s2_sparse localderivatives(numdofs);
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d0s2 localGrad;
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;

  // localDerivatives(i,j,k) = at node i, the derivative of j^{th} shape function w.r.t xi[k]
  getLocalDerivatives(&localderivatives, xi);

  for (int j = 0; j < numdofs/3; j++) {
    Node &nd = nodes[j]; 
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates; // dof contraction

  jacobian.getDeterminant(*jac);

  jacobian.getInverse(invjacobian);

  Tensor_d1s0 displacements(numdofs);
  disp.vectorToTensor(displacements);

  localGrad = localderivatives%displacements; // dof contraction

  //for(j = 0; j < numdofs; ++j)    
  //  fprintf(stderr,"disp[%d]=%e\n", j, displacements[j]);

  Tensor_d1s2_sparse * dgradUdqk = static_cast<Tensor_d1s2_sparse *>(_dgradUdqk);
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);

  (*gradU) = localGrad|invjacobian; // space contraction
  (*dgradUdqk) = localderivatives|invjacobian;

  //for (int q=0; q<9; q++)
  //  fprintf(stderr, "gradU[%d]=%e \n",q, (*gradU)[q]);

  //for (int q=0; q<8; q++)
  //  for (int p=0; p<9; p++)
  //    fprintf(stderr, "dgradUdqk[%d][%d]=%e \n",q, p, (*dgradUdqk)[q][p]);
}

void
ShapeFunction::getGradU(Tensor *_gradU, Node *nodes, double xi[3], Vector &disp)
{
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d1s2_sparse localderivatives(numdofs);

  getLocalDerivatives(&localderivatives, xi);

  for (int i = 0; i < numdofs; i++) {
    int j = i/3;
    Node &nd = nodes[j]; 
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates;//dof contraction

  jacobian.getInverse(invjacobian);

  Tensor_d0s2 localGrad;
  Tensor_d1s0 displacements(numdofs);
  // Can be done differently XFL
  disp.vectorToTensor(displacements);

  localGrad = (localderivatives%displacements);//dof contraction

  (*gradU) = (localGrad|invjacobian);//space contraction
}
