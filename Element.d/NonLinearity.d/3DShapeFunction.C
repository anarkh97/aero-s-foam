#include <stdio.h>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>


Tensor *
ShapeFunction::getGradUInstance(){

return (new Tensor_d0s2());
}



Tensor *
ShapeFunction::getDgradUDqkInstance(){

return (new Tensor_d1s2_sparse(numdofs));
}




void
ShapeFunction::getGlobalGrads(Tensor *_gradU, Tensor *_dgradUdqk, double *jac, Node *nodes, double xi[3],Vector &disp){





Tensor_d1s2_sparse localderivatives(numdofs);

Tensor_d1s0 nodescoordinates(numdofs);
Tensor_d0s2 localGrad;
Tensor_d0s2 jacobian;
Tensor_d0s2 invjacobian;

getLocalDerivatives(& localderivatives, xi);


for (int j = 0; j < numdofs/3; j++){
  Node &nd = nodes[j]; 
  nodescoordinates[3*j] = nd.x;
  nodescoordinates[3*j+1] = nd.y;
  nodescoordinates[3*j+2] = nd.z;
  
                              }


//isoparametric elements
jacobian = localderivatives%nodescoordinates;//dof contraction


jacobian.getDeterminant(*jac);
jacobian.getInverse(invjacobian);






Tensor_d1s0 displacements(numdofs);
disp.vectorToTensor(displacements);



localGrad = localderivatives%displacements;//dof contraction


//for(j = 0; j < numdofs; ++j)    
//    fprintf(stderr,"disp[%d]=%e\n", j, displacements[j]);

Tensor_d1s2_sparse * dgradUdqk = static_cast<Tensor_d1s2_sparse *>(_dgradUdqk);
Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);



(*gradU) = (localGrad|invjacobian);//space contraction
(*dgradUdqk) = (localderivatives|invjacobian);

//for (int q=0; q<9; q++)
//  fprintf(stderr, "gradU[%d]=%e \n",q, (*gradU)[q]);

//for (int q=0; q<8; q++)
// for (int p=0; p<9; p++)
//  fprintf(stderr, "dgradUdqk[%d][%d]=%e \n",q, p, (*dgradUdqk)[q][p]);
}




void
ShapeFunction::getGradU(Tensor *_gradU, Node *nodes, double xi[3],Vector &disp){

Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);
Tensor_d0s2 jacobian;
Tensor_d0s2 invjacobian;
Tensor_d1s0 nodescoordinates(numdofs);
Tensor_d1s2_sparse localderivatives(numdofs);

getLocalDerivatives(&localderivatives, xi);


for (int i = 0; i < numdofs; i++){
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






void
HexahedralShapeFunction::getLocalDerivatives(Tensor *_localDerivatives, double xi[3]){

Tensor_d1s2_sparse *localDerivatives  = static_cast<Tensor_d1s2_sparse *>(_localDerivatives);

static double localNodeCoordinates[24];
localNodeCoordinates[0]=-1;
localNodeCoordinates[1]=-1;
localNodeCoordinates[2]=-1;
localNodeCoordinates[3]=1;
localNodeCoordinates[4]=-1;
localNodeCoordinates[5]=-1;
localNodeCoordinates[6]=1;
localNodeCoordinates[7]=1;
localNodeCoordinates[8]=-1;
localNodeCoordinates[9]=-1;
localNodeCoordinates[10]=1;
localNodeCoordinates[11]=-1;
localNodeCoordinates[12]=-1;
localNodeCoordinates[13]=-1;
localNodeCoordinates[14]=1;
localNodeCoordinates[15]=1;
localNodeCoordinates[16]=-1;
localNodeCoordinates[17]=1;
localNodeCoordinates[18]=1;
localNodeCoordinates[19]=1;
localNodeCoordinates[20]=1;
localNodeCoordinates[21]=-1;
localNodeCoordinates[22]=1;
localNodeCoordinates[23]=1;





for (int k=0; k<8; k++) 
  for (int i=0; i<3; i++)  {         (*localDerivatives)[k][3*i]=(1.0/8.0)*(localNodeCoordinates[3*k]*(1+localNodeCoordinates[3*k+1]*xi[1])*(1+localNodeCoordinates[3*k+2]*xi[2]));
(*localDerivatives)[k][3*i+1]=(1/8.0)*(localNodeCoordinates[3*k+1]*(1+localNodeCoordinates[3*k]*xi[0])*(1+localNodeCoordinates[3*k+2]*xi[2]));
(*localDerivatives)[k][3*i+2]=(1/8.0)*(localNodeCoordinates[3*k+2]*(1+localNodeCoordinates[3*k]*xi[0])*(1+localNodeCoordinates[3*k+1]*xi[1]));
                           };

};
