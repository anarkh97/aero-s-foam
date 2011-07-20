#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/mathUtility.h>

FullSquareMatrix  
GaussIntgElement::stiffness(CoordSet& cs, double *k, int)
{
  int i;
  FullSquareMatrix kTan(numDofs(), k);

  Node *nodes = (Node *) dbg_alloca(numNodes()*sizeof(Node));
  int nnd = numNodes();
  int *ndn = (int *)dbg_alloca(nnd*sizeof(int)); 
  this->nodes(ndn); // PJSA 
  for(i = 0; i < nnd; ++i)
    nodes[i] = *cs[ndn[i]];
  int ndofs = numDofs();
  double *disp = (double *) dbg_alloca(ndofs*sizeof(double));
  for(i = 0; i < ndofs; ++i)
    disp[i] = 0;

  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradU = *shapeF->getGradUInstance();

  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqk = *shapeF->getDgradUDqkInstance();
 
  // NDofsx3x3x-> 6xNDofs
  Tensor &B = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs
  Tensor &DB = *strainEvaluator->getDBInstance(ndofs);

 
  Tensor &D = *strainEvaluator->getTMInstance();
  Tensor &e = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);
  
  //fprintf(stderr,"Je suis dans stiffness\n");

  int ngp = getNumGaussPoints();
  
  kTan.zero();
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);

    StackVector dispVec(disp,ndofs);
    shapeF->getGlobalGrads(&gradU, &dgradUdqk, &jac, nodes, point, dispVec);
    strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk);

    //material->getStress(&s,e,0);          
    //material->getTangentMaterial(&D, e, 0);
    //material->getStressAndTangentMaterial(&s, &D, e, 0);

    material->integrate(&s, &D, e, e, 0, 0, 0);

    //Tensor_d1s2_Ss23 & _B = static_cast<Tensor_d1s2_Ss23 &>(B);
    //for (int k = 0; k<24; ++k)
    // for (int j = 0; j<6; ++j) 
    //    fprintf(stderr,"B[%d][%d]=%e\n", k, j, _B[k][j]);
              
    temp1 =  D || B;
    temp2 =   B||temp1;
    temp3 = DB || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jac))*temp3;

    kTan += temp3;
  }
  delete &temp1;
  delete &s;
  delete &gradU;
  delete &dgradUdqk;
  delete &B;
  delete &DB; 
  delete &e;
  delete &D;
  //fprintf(stderr, "Diag:\n");
  //for(i = 0; i< ndofs; ++i)
  //    fprintf(stderr, "%e ", kTan[i][i]);
  //fprintf(stderr, "\n");
  return kTan;
}

FullSquareMatrix  
GaussIntgElement::massMatrix(CoordSet& cs, double* k, int)
{
  FullSquareMatrix m(numDofs(), k);
  cerr << "GaussIntgElement::massMatrix not implemented\n"; exit(-1);
  return m;
}
 
void
GaussIntgElement::getStiffAndForce(Node *nodes, double *disp,
                                   double *state, FullSquareMatrix &kTan,
                                   double *force)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradU = *shapeF->getGradUInstance();

  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqk = *shapeF->getDgradUDqkInstance();
  
  // NDofsx3x3x-> 6xNDofs
  Tensor &B = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but diagonal in dofs
  Tensor &DB = *strainEvaluator->getDBInstance(ndofs);

  Tensor &D = *strainEvaluator->getTMInstance();
  Tensor &e = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  Tensor_d1s0 nodeforce(ndofs);
  Tensor_d1s0 temp0(ndofs);
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);

  fprintf(stderr,"Je suis dans getStiffAndForce\n");

  int i,j;
  int ngp = getNumGaussPoints();
  
  kTan.zero();
  
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);

    StackVector dispVec(disp,ndofs);
        
    shapeF->getGlobalGrads(&gradU, &dgradUdqk,  &jac, nodes, point, dispVec);
    strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk);

    //material->getStress(&s, e, 0);        
    //material->getStressAndTangentMaterial(&s, &D, e, 0);
        
    material->integrate(&s, &D,  e, e,0, 0, 0);

    //Tensor_d0s4_Ss12s34 & _D = static_cast<Tensor_d0s4_Ss12s34 &>(D);
    //for (int k = 0; k<6; ++k)
    //  for (int j = 0; j<6; ++j) 
    //    fprintf(stderr,"D[%d][%d]=%e\n", k, j, _D[k][j]);

    //Tensor_d1s2_Ss23 & _B = static_cast<Tensor_d1s2_Ss23 &>(B);
    //for (int k = 0; k<24; ++k)
    //  for (int j = 0; j<6; ++j) 
    //    fprintf(stderr,"B[%d][%d]=%e\n", k, j, _B[k][j]);

    temp0 = s || B;
    temp0 = (weight*jac)*temp0;
    nodeforce = nodeforce + temp0;
    temp1 =  D || B;
    temp2 =   B||temp1;
    temp3 = DB || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jac))*temp3;
    kTan += temp3;
    temp3 = DB || s;

    //Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(s);
    //for (int ii=0; ii<6; ++ii)
    //  fprintf(stderr," %e", stress[ii]);
    //fprintf(stderr, "\n");
    /*   fprintf(stderr, "D || B\n");
    Tensor_d2s0 &_temp3 = static_cast<Tensor_d2s0 &>(temp3);
    for (int k = 0; k<24; ++k) {
      for (int j = 0; j<24; ++j) 
        fprintf(stderr," %e", _temp3[k*24+j]);
      fprintf(stderr, "\n");
     }*/
  }

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];
    //fprintf(stderr,"force[%d]=%e\n", j, force[j]);
  }
  // fprintf(stderr, "kTan\n");
  // kTan.print();

  ////////////////////////////////////////////////////////
  //Check convergence with approximate tangent stiffness//
  ////////////////////////////////////////////////////////

  Tensor &eleft = *strainEvaluator->getStrainInstance();
  Tensor &eright = *strainEvaluator->getStrainInstance();
  Tensor &DBleft = *strainEvaluator->getDBInstance(ndofs);
  Tensor &DBright = *strainEvaluator->getDBInstance(ndofs);
  Tensor &gradUleft = *shapeF->getGradUInstance();
  Tensor &gradUright = *shapeF->getGradUInstance();
  Tensor &dgradUdqkleft = *shapeF->getDgradUDqkInstance();
  Tensor &dgradUdqkright = *shapeF->getDgradUDqkInstance();
  Tensor &sleft = *strainEvaluator->getStressInstance();
  Tensor &sright = *strainEvaluator->getStressInstance();
  Tensor &Bleft = *strainEvaluator->getBInstance(ndofs);
  Tensor &Bright = *strainEvaluator->getBInstance(ndofs);
  
  double *kt = new double[ndofs*ndofs];
  FullSquareMatrix kTant(ndofs, kt);
 
  int a,b;
  double dispsqNorm =0.0;
  for (b=0; b < ndofs; ++b)
    dispsqNorm+=disp[b]*disp[b];
  double epsilon = 1e-8; //(1e-8)*sqrt(dispsqNorm);
  //if(dispsqNorm == 0.0) epsilon = 1e-8;
  kTant.zero();
  fprintf(stderr, "epsilon = %e\n" , epsilon);
  double* dispLeft = new double[ndofs];
  double* dispRight = new double[ndofs];
  for (b=0; b < ndofs; ++b) {
    dispLeft[b] = disp[b];
    dispRight[b] = disp[b]; 
  }

  for (b=0; b < ndofs; ++b) {

    if(disp[b] > 1 || disp[b] < -1)
      epsilon = 1e-8*disp[b];
    else
      epsilon = 1e-8;
    dispLeft[b] = disp[b]+epsilon;
    dispRight[b] = disp[b]-epsilon;
    Tensor_d1s0 nodeforceLeft(ndofs);
    Tensor_d1s0 nodeforceRight(ndofs);
    StackVector dispVecTestLeft(dispLeft,ndofs);
    StackVector dispVecTestRight(dispRight,ndofs);
     
    for(i = 0; i < ngp; i++) {
      double point[3], weight, jac;
      getGaussPointAndWeight(i, point, weight);

      shapeF->getGlobalGrads(&gradUleft, &dgradUdqkleft,  &jac, nodes, point, dispVecTestLeft);
      strainEvaluator->getEBandDB(eleft, Bleft, DBleft, gradUleft, dgradUdqkleft);
      material->integrate(&sleft, &D,  eleft, eleft,0, 0, 0);
      temp0 = sleft || Bleft;
      temp0 = (weight*jac)*temp0;
      nodeforceLeft = nodeforceLeft + temp0;

      shapeF->getGlobalGrads(&gradUright, &dgradUdqkright,  &jac, nodes, point, dispVecTestRight);
      strainEvaluator->getEBandDB(eright, Bright, DBright, gradUright, dgradUdqkright);
      material->integrate(&sright, &D,  eright, eright,0, 0, 0);
      temp0 = sright || Bright;
      temp0 = (weight*jac)*temp0;
      nodeforceRight = nodeforceRight + temp0;
    }

    for (a=0; a < ndofs; ++a) {              
      kTant[a][b] = (nodeforceLeft[a] - nodeforceRight[a])/(2*epsilon);
      //fprintf(stderr, "kTant[%d][%d] = %e\n" ,a ,b , kTant[a][b]);                                  
    }

    /* if(b == 1) fprintf(stderr, "%e vs %e or %e or %e\n", kTan[b][b], kTant[b][b],
                          (-nodeforceLeft[b]-force[b])/epsilon, (force[b]+nodeforceRight[b])/epsilon);*/

    dispLeft[b] = disp[b];
    dispRight[b] = disp[b];              
  }
  //double diffcol[ndofs];
  //double col[ndofs];
  double diffsqNorm = 0.0;
  double sqNorm = 0.0;
  double relativeNorm;

  for (a=0; a < ndofs; ++a)         
    for (b=0; b < ndofs; ++b) {
      // diffcol[a]+=fabs(kTan[a][b]-kTant[a][b]);
      //col[a]+=fabs(kTant[a][b])                   

      diffsqNorm+=(kTan[a][b]-kTant[a][b])*(kTan[a][b]-kTant[a][b]);
      sqNorm+=kTant[a][b]*kTant[a][b];
    }
  relativeNorm = sqrt(diffsqNorm/sqNorm);
  fprintf(stderr, "Relative Norm = %e\n" , relativeNorm );
  //fprintf(stderr, "%e vs %e\n", kTan[0][0], kTant[0][0]);
        
  delete &Bleft;
  delete &Bright;        
  delete &sleft;
  delete &sright;
  delete &gradUleft;
  delete &gradUright;
  delete &dgradUdqkleft; 
  delete &dgradUdqkright;
  delete &eleft;
  delete &eright;
  delete &DBleft;
  delete &DBright;

  /////////////////////////////////////////////////////////////////////////////

  delete &temp1;
  delete &gradU;
  delete &dgradUdqk;
  delete &B;
  delete &DB;
  delete &e;
  delete &s;
  delete &D;
}

//void
//GaussIntgElement::updateStates(Node *nodes, double *state, double *un,double *unp){}
/*
void
GaussIntgElement::updateStates(Node *nodes, double *state, double *un,double *unp){

StrainEvaluator *strainEvaluator = getStrainEvaluator();
ShapeFunction *shapeF = getShapeFunction();
NLMaterial *material = getMaterial();

Tensor &gradUn = *shapeF->getGradUInstance();
Tensor &gradUnp = *shapeF->getGradUInstance();


Tensor &en = *strainEvaluator->getStrainInstance();
Tensor &enp = *strainEvaluator->getStrainInstance();


int ndofs = numDofs();
int ngp  = getNumGaussPoints();
int nstatepgp = material->getNumStates();

double point[3], weight;

StackVector dispn(un,ndofs);
StackVector dispnp(unp,ndofs);


for (int i=0; i<ngp; ++i){ 

  getGaussPointAndWeight(i, point, weight);
  shapeF->getGradU(&gradUn, nodes, point, dispn);
  shapeF->getGradU(&gradUnp, nodes, point, dispnp);
  strainEvaluator->getE(en, gradUn);
  strainEvaluator->getE(enp, gradUnp);
  material->updateStates(en, enp, state + nstatepgp*i);
                         };
};
*/


void 
GaussIntgElement::integrate(Node *nodes, double *dispn,  double *staten,
                            double *dispnp, double *statenp,
                            FullSquareMatrix &kTan,
                            double *force, double)
{
  //for(int i=0; i<8; ++i) cerr << nodes[i].x+dispn[3*i+0] << "," << nodes[i].y+dispn[3*i+1] << "," << nodes[i].z+dispn[3*i+2] << " ";
  //cerr << endl;
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();
  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqkn = *shapeF->getDgradUDqkInstance();
  Tensor &dgradUdqknp = *shapeF->getDgradUDqkInstance();

  // NDofsx3x3x-> 6xNDofs
  Tensor &Bn = *strainEvaluator->getBInstance(ndofs);
  Tensor &Bnp = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but sparse
  Tensor &DBn = *strainEvaluator->getDBInstance(ndofs);
  Tensor &DBnp = *strainEvaluator->getDBInstance(ndofs);

  Tensor &Dnp = *strainEvaluator->getTMInstance();
  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();
  
  Tensor_d1s0 nodeforce(ndofs);
  Tensor_d1s0 temp0(ndofs);
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs, false);
  Tensor_d2s0 temp3(ndofs, false);

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();
  
  kTan.zero();

  //fprintf(stderr,"Je suis dans integrate\n");
  
  for(i = 0; i < ngp; i++) {

    double point[3], weight, jacn, jacnp;
    StackVector dispVecn(dispn, ndofs);
    StackVector dispVecnp(dispnp, ndofs); 
 
    getGaussPointAndWeight(i, point, weight);
/* PJSA seems to be unnecessary
    shapeF->getGradU(&gradUn, nodes, point, dispVecn);
    shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
*/
/* PJSA seems to be unnecessary
    strainEvaluator->getE(en, gradUn);
    strainEvaluator->getE(enp, gradUnp);  
*/
    shapeF->getGlobalGrads(&gradUn, &dgradUdqkn, &jacn, nodes, point, dispVecn);
    shapeF->getGlobalGrads(&gradUnp, &dgradUdqknp, &jacnp, nodes, point, dispVecnp);

    strainEvaluator->getEBandDB(en, Bn, DBn, gradUn, dgradUdqkn);
    strainEvaluator->getEBandDB(enp, Bnp, DBnp, gradUnp, dgradUdqknp);

    //material->updateStates(en, enp, state + nstatepgp*i);
    //material->getStress(&s, e, 0);       
    //material->getStressAndTangentMaterial(&s, &D, enp, 0);
    material->integrate(&s, &Dnp, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, 0);

    //std::cerr << "s = "; s.print();
    //std::cerr << "Dnp = "; Dnp.print();
    temp0 = s || Bnp;
    temp0 = (weight*jacnp)*temp0;
    nodeforce = nodeforce + temp0;
    temp1 =  Dnp || Bnp;
    temp2 =   Bnp||temp1;
    //cerr << "DBnp = "; DBnp.print();
    temp3 = DBnp || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jacnp))*temp3;
    kTan += temp3;

/*
    temp1 =  D || B;
    temp2 =   B||temp1;
    temp3 = DB || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jac))*temp3;

    kTan += temp3;
*/
  }

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];}

  delete &temp1;
  delete &gradUn;
  delete &gradUnp;
  delete &dgradUdqkn;
  delete &dgradUdqknp;
  delete &Bn;
  delete &DBn;
  delete &Bnp;
  delete &DBnp;
  delete &en;
  delete &enp;
  delete &s;
  delete &Dnp;

  //cerr << "K = "; kTan.print();
}

void
GaussIntgElement::initStates(double *st)
{
  NLMaterial *material = getMaterial();
  int ninterns = material->getNumStates();
  int ngp = getNumGaussPoints();
 
  for(int i = 0; i < ngp; ++i)
    material->initStates(st+i*ninterns);
}

void
GaussIntgElement::updateStates(Node *nodes, double *state, double *dispn, double *dispnp)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();
  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqkn = *shapeF->getDgradUDqkInstance();
  Tensor &dgradUdqknp = *shapeF->getDgradUDqkInstance();

  // NDofsx3x3x-> 6xNDofs
  Tensor &Bn = *strainEvaluator->getBInstance(ndofs);
  Tensor &Bnp = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but sparse
  Tensor &DBn = *strainEvaluator->getDBInstance(ndofs);
  Tensor &DBnp = *strainEvaluator->getDBInstance(ndofs);

  Tensor &Dnp = *strainEvaluator->getTMInstance();
  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  int nstatepgp = material->getNumStates();

  for(int i = 0; i < getNumGaussPoints(); i++) {

    double point[3], weight, jacn, jacnp;
    StackVector dispVecn(dispn, ndofs);
    StackVector dispVecnp(dispnp, ndofs);

    getGaussPointAndWeight(i, point, weight);

    //shapeF->getGradU(&gradUn, nodes, point, dispVecn);
    //shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
    //strainEvaluator->getE(en, gradUn);
    //strainEvaluator->getE(enp, gradUnp);  
    shapeF->getGlobalGrads(&gradUn, &dgradUdqkn, &jacn, nodes, point, dispVecn);
    shapeF->getGlobalGrads(&gradUnp, &dgradUdqknp, &jacnp, nodes, point, dispVecnp);
    strainEvaluator->getEBandDB(en, Bn, DBn, gradUn, dgradUdqkn);
    strainEvaluator->getEBandDB(enp, Bnp, DBnp, gradUnp, dgradUdqknp);

    //material->updateStates(en, enp, state + nstatepgp*i);
    //material->getStress(&s, e, 0);       
    //material->getStressAndTangentMaterial(&s, &D, enp, 0);
    double *state_copy = new double[nstatepgp];
    for(int j = 0; j < nstatepgp; ++j) state_copy[j] = state[nstatepgp*i+j];
    material->integrate(&s, &Dnp, en, enp,
                        state_copy, state + nstatepgp*i, 0);
    delete [] state_copy;
  }

  delete &gradUn;
  delete &gradUnp;
  delete &dgradUdqkn;
  delete &dgradUdqknp;
  delete &Bn;
  delete &DBn;
  delete &Bnp;
  delete &DBnp;
  delete &en;
  delete &enp;
  delete &s;
  delete &Dnp;
}

