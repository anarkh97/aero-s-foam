#ifndef _GAUSSINTGELEM_H_
#define _GAUSSINTGELEM_H_
#include <alloca.h>
#include <Math.d/FullSquareMatrix.h>
#include <Element.d/NLElement.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Utils.d/NodeSpaceArray.h>
#include <Math.d/TTensor.h>
#include <Utils.d/dbg_alloca.h>

class StrainEvaluator;
template <class TT> class GenStrainEvaluator;

class GaussIntgElement : public MatNLElement
{
  protected:
    virtual int getNumGaussPoints() = 0;
    virtual void getGaussPointAndWeight(int, double *, double &) = 0;
    virtual void getLocalNodalCoords(int, double *) = 0;
    virtual ShapeFunction *getShapeFunction() = 0;
    virtual StrainEvaluator *getStrainEvaluator() = 0;
    virtual NLMaterial *getMaterial() = 0;

  public:
    void getStiffAndForce(Node *nodes, double *disp,
                          double *state, FullSquareMatrix &kTan,
                          double *force);
    FullSquareMatrix  stiffness(CoordSet& cs, double *k, int flg=1);
    FullSquareMatrix massMatrix(CoordSet& cs, double *m, int flg=1);
    void updateStates(Node *node, double *state, double *un, double *unp);
    template <class MatrixType, class MaterialType>
      void integrate(const MaterialType &,
                     Node *nodes, double *dispn, double *staten,
                     double *dispnp, double *statenp,
                     MatrixType &kTan, double *force, double dt=0.0);
    void integrate(Node *nodes, double *dispn, double *staten,
                   double *dispnp, double *statenp,
                   FullSquareMatrix &kTan, double *force, double dt=0.0);
    int numStates() {
      int nGP = getNumGaussPoints();
      NLMaterial *mat = getMaterial();
      int nsGP = mat->getNumStates(); 
      return nGP*nsGP;
    }
    void initStates(double *);
    // the following functions return result for postprocessing at every node
    void getStrainTens(Node *nodes, double *dispnp, double (*result)[9]);
    void getVonMisesStrain(Node *nodes, double *dispnp, double *result);
    void getStressTens(Node *nodes, double *dispn, double *staten,
                       double *dispnp, double *statenp, double (*result)[9]);
    void getVonMisesStress(Node *nodes, double *dispn, double *staten,
                           double *dispnp, double *statenp, double *result);
    void getEquivPlasticStrain(double *statenp, double *result);
};

template <class TensorTypes>
class GenGaussIntgElement : public MatNLElement
{
  protected:
    virtual int getNumGaussPoints() = 0;
    virtual void getGaussPointAndWeight(int, double *, double &) = 0;
    virtual GenShapeFunction<TensorTypes> *getShapeFunction() = 0;
    virtual GenStrainEvaluator<TensorTypes> *getGenStrainEvaluator() = 0;
    virtual NLMaterial *getMaterial() = 0;

   public:
    void getStiffAndForce(Node *nodes, double *disp,
                          double *state, FullSquareMatrix &kTan,
                          double *force);
    FullSquareMatrix  stiffness(CoordSet& cs, double *k, int flg=1);
    FullSquareMatrix massMatrix(CoordSet& cs, double *m, int flg=1);
    void updateStates(Node *node, double *state, double *un, double *unp);
    template <class MatrixType, class MaterialType>
      void integrate(const MaterialType &,
                     Node *nodes, double *dispn, double *staten,
                     double *dispnp, double *statenp,
                     MatrixType &kTan, double *force, double dt=0.0);
    void integrate(Node *nodes, double *dispn, double *staten,
                   double *dispnp, double *statenp,
                   FullSquareMatrix &kTan, double *force, double dt=0.0);
    int numStates() {
      int ngp = getNumGaussPoints();
      NLMaterial *mat = getMaterial();
      int nst = mat->getNumStates();
      return nst*ngp;
    }
    void initStates(double *);
};

template <class TensorTypes>
FullSquareMatrix  
GenGaussIntgElement<TensorTypes>::stiffness(CoordSet& cs, double *k, int)
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
  double *disp = (double *)dbg_alloca(ndofs*sizeof(double));
  for(i = 0; i < ndofs; ++i)
    disp[i] = 0;

  GenShapeFunction<TensorTypes> *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  GenStrainEvaluator<TensorTypes> *strainEvaluator = getGenStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

//  fprintf(stderr, "Material is %p\n", material);

  // Obtain the storage for gradU ( 3xndim )
  typename TensorTypes::GradUTensor gradU;

  // Obtain the storage for dgradUdqk ( ndof x3x ndim )
  typename TensorTypes::GradUDerivTensor dgradUdqk;
 
  // NDofsx3x3x-> 6xNDofs
  typename TensorTypes::BTensor B;

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs
  typename TensorTypes::DBTensor DB;

  typename TensorTypes::DTensor D;
  typename TensorTypes::StressTensor e;
  typename TensorTypes::StressTensor s;
  typename TensorTypes::BTensor temp1;
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);
  
  //fprintf(stderr,"Je suis dans stiffness\n");

  int ngp = getNumGaussPoints();
 
  kTan.zero();
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);
    // fprintf(stderr, "Material %p\n", material);

    StackVector dispVec(disp,ndofs);
    // fprintf(stderr, "Material 2 %p\n", material);
    shapeF->getGlobalGrads(gradU, dgradUdqk, &jac, nodes, point, dispVec);
    // fprintf(stderr, "Material 3 %p\n", material);
    //fprintf(stderr, "dgradUdqk:\n");
    //int kkk;
    //for(kkk = 0; kkk < 9; ++kkk)
    //  fprintf(stderr,"%e %e %e    %e %e %e\n", dgradUdqk[kkk][0][0], dgradUdqk[kkk][0][1],
    //          dgradUdqk[kkk][0][2],
    //          dgradUdqk[kkk][1][0], dgradUdqk[kkk][1][1], dgradUdqk[kkk][1][2]);
    strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk);

    //material->getStress(&s,e,0);          
    //material->getTangentMaterial(&D, e, 0);
    //material->getStressAndTangentMaterial(&s, &D, e, 0);

    //fprintf(stderr, "About to integrate! %p\n", material);
    material->integrate(&s, &D, e, e, 0, 0, 0);
    //fprintf(stderr, "Have integrated\n");

    //Tensor_d1s2_Ss23 & _B = static_cast<Tensor_d1s2_Ss23 &>(B);
    //for (int k = 0; k<24; ++k)
    // for (int j = 0; j<6; ++j) 
    //    fprintf(stderr,"B[%d][%d]=%e\n", k, j, _B[k][j]);
              
    // temp1 =  D || B;
    temp1 = dblContractTransp(D,B);
    //fprintf(stderr, "D: %e %e %e %e %e %e %e %e %e\n", D[0][0], D[0][1], D[0][2], D[1][0], D[1][1], D[1][2],
    //          D[2][0], D[2][1], D[2][2]);
    //fprintf(stderr, "dSdqk:\n");
    //for(kkk = 0; kkk < 9; ++kkk)
    //  fprintf(stderr,"%e %e %e\n", temp1[kkk][0], temp1[kkk][1], temp1[kkk][2]);
    //fprintf(stderr, "B:\n");
    //for(kkk=0; kkk < ndofs;++kkk)
    //  fprintf(stderr, "%d %e %e %e\n", kkk, B[kkk][0], B[kkk][1], B[kkk][2]);
    temp2 =   B||temp1;
    // We should test if the strain is non linear to call the following two lines
    if(strainEvaluator->isNonLinear()) {
      SimpleTensor<SimpleTensor<double,TensorTypes::ndofs>,TensorTypes::ndofs> kg;
      kg=DB || s;
      for(int j=0; j < TensorTypes::ndofs; ++j)
        for(int k=0; k < TensorTypes::ndofs; ++k)
          temp3[j*TensorTypes::ndofs+k] = kg[j][k];
      temp2 = temp3 + temp2;
    }
    //fprintf(stderr, "Weight %e jac %e\n",weight,jac);
    //for(kkk = 0; kkk < 9; ++kkk) {
    //  for(int i=0; i < 9; ++i)
    //    fprintf(stderr, "%e ", temp2[kkk*9+i]);
    //  fprintf(stderr, "\n");
    //}
    temp3 = (weight * fabs(jac))*temp2;
    //fprintf(stderr, "K\n");
    //for(kkk = 0; kkk < 9; ++kkk) {
    //  for(int i=0; i < 9; ++i)
    //    fprintf(stderr, "%e ", temp3[kkk*9+i]);
    //  fprintf(stderr, "\n");
    //}

    kTan += temp3;
  }
  //fprintf(stderr, "Diag:\n");
  //for(i = 0; i< ndofs; ++i)
  //    fprintf(stderr, "%e ", kTan[i][i]);
  //fprintf(stderr, "\n");
  return kTan;
}

template <class TensorType>
FullSquareMatrix  
GenGaussIntgElement<TensorType>::massMatrix(CoordSet&, double *m, int)
{
  return FullSquareMatrix(numDofs(), m);
}
 
template <class TensorType>
void
GenGaussIntgElement<TensorType>::getStiffAndForce(Node *nodes, double *disp,
                                                  double *state, FullSquareMatrix &kTan,
                                                  double *force)
{
  
  int ndofs = numDofs();
  // ShapeFunction *shapeF = getShapeFunction();
  GenShapeFunction<TensorType> *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  //Tensor &gradU = *shapeF->getGradUInstance();
  typename TensorType::GradUTensor gradU;

  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  typename TensorType::GradUDerivTensor dgradUdqk;
  
  // NDofsx3x3x-> 6xNDofs
  typename TensorType::BTensor B;

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but diagonal in dofs
  typename TensorType::DBTensor DB;

  typename TensorType::DTensor D;
  typename TensorType::StrainTensor e;
  typename TensorType::StressTensor s;

  SimpleTensor<double, TensorType::ndofs> temp0, nodeforce;
  typename TensorType::BTensor temp1;
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);

  //fprintf(stderr,"Je suis dans getStiffAndForce\n");

  int i,j;
  int ngp = getNumGaussPoints();
  
  kTan.zero();
  
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);

    StackVector dispVec(disp,ndofs);
        
    shapeF->getGlobalGrads(gradU, dgradUdqk,  &jac, nodes, point, dispVec);
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
    temp1 =  dblContractTransp(D, B);
    temp2 =  B||temp1;
    if(strainEvaluator->isNonLinear()) {
      SimpleTensor<SimpleTensor<double,TensorType::ndofs>,TensorType::ndofs> kg;
      kg=DB || s;
      for(int j=0; j < TensorType::ndofs; ++j)
        for(int k=0; k < TensorType::ndofs; ++k)
          temp3[j*TensorType::ndofs+k] = kg[j][k];
        temp3 = temp3 + temp2;
      }
      temp3 = (weight * fabs(jac))*temp3;
      kTan += temp3;

      //Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(s);
      //for (int ii=0; ii<6; ++ii)
      //  fprintf(stderr," %e", stress[ii]);
      //fprintf(stderr, "\n");
      /*fprintf(stderr, "D || B\n");
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

  typename TensorType::StrainTensor eleft, eright, sleft, sright;
  typename TensorType::DBTensor DBleft, DBright;
  typename TensorType::BTensor Bleft, Bright;
  typename TensorType::GradUTensor gradUleft;
  typename TensorType::GradUTensor gradUright;
  typename TensorType::GradUDerivTensor dgradUdqkleft;
  typename TensorType::GradUDerivTensor dgradUdqkright;
  
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
  for (b=0; b < ndofs; ++b){
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
    SimpleTensor<double, TensorType::ndofs> nodeforceLeft, nodeforceRight;
    StackVector dispVecTestLeft(dispLeft,ndofs);
    StackVector dispVecTestRight(dispRight,ndofs);
      
    for(i = 0; i < ngp; i++) {
      double point[3], weight, jac;
      getGaussPointAndWeight(i, point, weight);

      shapeF->getGlobalGrads(gradUleft, dgradUdqkleft,  &jac, nodes, point, dispVecTestLeft);
      strainEvaluator->getEBandDB(eleft, Bleft, DBleft, gradUleft, dgradUdqkleft);
      material->integrate(&sleft, &D,  eleft, eleft,0, 0, 0);
      temp0 = sleft || Bleft;
      temp0 = (weight*jac)*temp0;
      nodeforceLeft = nodeforceLeft + temp0;

      shapeF->getGlobalGrads(gradUright, dgradUdqkright,  &jac, nodes, point, dispVecTestRight);
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
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::updateStates(Node *nodes, double *state, double *un, double *unp){}
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


template <class TensorType>
void 
GenGaussIntgElement<TensorType>::integrate(Node *nodes, double *dispn,  double *staten,
                                           double *dispnp, double *statenp,
                                           FullSquareMatrix &kTan,
                                           double *force, double)
{
  int ndofs = numDofs();
  GenShapeFunction<TensorType> *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  typename TensorType::GradUTensor gradUn;
  typename TensorType::GradUTensor gradUnp;
  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  typename TensorType::GradUDerivTensor dgradUdqkn;
  typename TensorType::GradUDerivTensor dgradUdqknp;

  // NDofsx3x3x-> 6xNDofs
  typename TensorType::BTensor Bn;
  typename TensorType::BTensor Bnp;

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but sparse
  typename TensorType::DBTensor DBn;
  typename TensorType::DBTensor DBnp;

  typename TensorType::DTensor Dnp;
  typename TensorType::StressTensor en;
  typename TensorType::StressTensor enp;
  typename TensorType::StressTensor s;
  
  SimpleTensor<double,TensorType::ndofs> temp0;
  SimpleTensor<double,TensorType::ndofs> nodeforce;
  typename TensorType::BTensor temp1;
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();
  
  kTan.zero();
  nodeforce = 0;
  
  for(i = 0; i < ngp; i++) {

    double point[3], weight, jacn, jacnp;
    StackVector dispVecn(dispn,ndofs);
    StackVector dispVecnp(dispnp,ndofs); 

    getGaussPointAndWeight(i, point, weight);

    /*shapeF->getGradU(gradUn, nodes, point, dispVecn);
    shapeF->getGradU(gradUnp, nodes, point, dispVecnp);

    strainEvaluator->getE(en, gradUn);
    strainEvaluator->getE(enp, gradUnp);  */

    shapeF->getGlobalGrads(gradUn, dgradUdqkn,  &jacn, nodes, point, dispVecn);
    shapeF->getGlobalGrads(gradUnp, dgradUdqknp,  &jacnp, nodes, point, dispVecnp);

    strainEvaluator->getEBandDB(en, Bn, DBn, gradUn, dgradUdqkn);
    strainEvaluator->getEBandDB(enp, Bnp, DBnp, gradUnp, dgradUdqknp);

    //material->updateStates(en, enp, state + nstatepgp*i);
    //material->getStress(&s, e, 0);       
    //material->getStressAndTangentMaterial(&s, &D, enp, 0);

    material->integrate(&s, &Dnp, en, enp,
                        staten + nstatepgp*i,statenp + nstatepgp*i , 0);

    temp0 = s || Bnp;
    temp0 = (weight*jacnp)*temp0;
    nodeforce = nodeforce + temp0;
    //temp1 =  Dnp || Bnp;
    temp1 = dblContractTransp(Dnp, Bnp);
    temp2 =   Bnp||temp1;
    if(strainEvaluator->isNonLinear()) {
      SimpleTensor<SimpleTensor<double,TensorType::ndofs>,TensorType::ndofs> kg;
      kg=DBnp || s;
      for(int j=0; j < TensorType::ndofs; ++j)
        for(int k=0; k < TensorType::ndofs; ++k)
          temp3[j*TensorType::ndofs+k] = kg[j][k];
        temp2 = temp3 + temp2;
    }
    temp3 = (weight * fabs(jacnp))*temp2;
    kTan += temp3;
  }

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];}
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::initStates(double *st)
{
  NLMaterial *material = getMaterial();
  int ninterns = material->getNumStates();
  int ngp = getNumGaussPoints();
 
  for(int i = 0; i < ngp; ++i)
    material->initStates(st+i*ninterns);
}

#endif
