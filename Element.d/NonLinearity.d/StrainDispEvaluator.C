#include <Math.d/matrix.h> 
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainDispEvaluator.h>
#include <stdio.h>


LinearStrain linearStrain;
GreenLagrangeStrain greenLagrangeStrain;


Tensor *
LinearStrain::getTMInstance(){
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}


Tensor *
LinearStrain::getStressInstance(){
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}


Tensor *
LinearStrain::getStrainInstance(){
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}


Tensor *
LinearStrain::getBInstance(int numdofs){
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}


Tensor *
LinearStrain::getDBInstance(int numdofs){
  Tensor_d2s2_Sd12s34_null *DB = new Tensor_d2s2_Sd12s34_null(numdofs);
  return DB;
}





void 
LinearStrain::getEBandDB(Tensor &_e, Tensor & bB, Tensor &DB,
                         Tensor &_gradU,Tensor &_dgradUdqk){

Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);
Tensor_d1s2_sparse & dgradUdqk = static_cast<Tensor_d1s2_sparse &>(_dgradUdqk);
Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(bB);
Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);


//int size = B.getSize();

Tensor_d0s2 enonsym;
Tensor_d0s2 tgradU;


gradU.getTranspose(tgradU);


enonsym = (1/2.)*(gradU + tgradU);

enonsym.convertToSym(e);



B = dgradUdqk.symPart();

};


void 
LinearStrain::getE(Tensor &_e,Tensor &_gradU){

Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);
Tensor_d0s2 enonsym;
Tensor_d0s2 tgradU;

gradU.getTranspose(tgradU);

enonsym = (1/2.)*(gradU + tgradU);
enonsym.convertToSym(e);



};




Tensor *
GreenLagrangeStrain::getTMInstance(){
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}



Tensor *
GreenLagrangeStrain::getStressInstance(){
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}


Tensor *
GreenLagrangeStrain::getStrainInstance(){
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
GreenLagrangeStrain::getBInstance(int numdofs){
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}



Tensor *
GreenLagrangeStrain::getDBInstance(int numdofs){
  Tensor_d2s2_Sd12s34_sparse *DB = new Tensor_d2s2_Sd12s34_sparse(numdofs);
  return DB;
}



void 
GreenLagrangeStrain::getEBandDB(Tensor &_e, Tensor &bB, Tensor &_DB,Tensor &_gradU, Tensor &_dgradUdqk){

Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);
Tensor_d1s2_sparse & dgradUdqk = static_cast<Tensor_d1s2_sparse &>(_dgradUdqk);
Tensor_d2s2_Sd12s34_sparse & DB = static_cast<Tensor_d2s2_Sd12s34_sparse &>(_DB);
Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(bB);
Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

int size = B.getSize();
Tensor_d0s2 enonsym;
Tensor_d0s2 tgradU;

gradU.getTranspose(tgradU);

Tensor_d0s2 temp1;


temp1 = tgradU|gradU;


enonsym = (1/2.)*(tgradU + (gradU + temp1));

enonsym.convertToSym(e);

Tensor_d1s2_full temp2(size);
temp2 = tgradU | dgradUdqk;

B = dgradUdqk.symPart() + temp2.symPart();


dgradUdqk.getSymSquare(DB);

};



void 
GreenLagrangeStrain::getE(Tensor &_e,Tensor &_gradU){

Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

Tensor_d0s2 enonsym;
Tensor_d0s2 tgradU;

gradU.getTranspose(tgradU);

Tensor_d0s2 temp;
temp = tgradU|gradU;

enonsym = (1/2.)*(tgradU + (gradU + temp));
enonsym.convertToSym(e);


};
