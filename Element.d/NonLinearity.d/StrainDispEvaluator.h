#ifndef _STRAINDISPEVALUATOR_H_
#define _STRAINDISPEVALUATOR_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>



//
//Computes the geometrical part of the stiffness matrix 


class StrainEvaluator {

public:
        virtual Tensor *getTMInstance() = 0;
	virtual Tensor *getStressInstance() = 0;
        virtual Tensor *getStrainInstance() = 0;
        virtual Tensor *getBInstance(int numdofs) = 0;
	virtual Tensor *getDBInstance(int numdofs) = 0;
        virtual void getEBandDB(Tensor &e, Tensor &B, Tensor &DB,Tensor &gradU, Tensor &dgradUdqk) = 0;
        virtual void getE(Tensor &e, Tensor &gradU) = 0;

};



class LinearStrain : public StrainEvaluator {

public:
	Tensor *getTMInstance();
        Tensor *getStressInstance();
        Tensor *getStrainInstance();
        Tensor *getBInstance(int numdofs);
	Tensor *getDBInstance(int numdofs);
        void getEBandDB(Tensor &e, Tensor &B, Tensor &DB,Tensor &gradU, Tensor &dgradUdqk);
        void getE(Tensor &e, Tensor &gradU);




};



class GreenLagrangeStrain : public StrainEvaluator {

public:
        Tensor *getTMInstance();
	Tensor *getStressInstance();
        Tensor *getStrainInstance();
        Tensor *getBInstance(int numdofs);
	Tensor *getDBInstance(int numdofs);
        void getEBandDB(Tensor &e, Tensor &B, Tensor &DB,Tensor &gradU, Tensor &dgradUdqk);
        void getE(Tensor &e, Tensor &gradU);
};

template<class TensorTypes>
class GenStrainEvaluator {
 public:
  virtual void getEBandDB(typename TensorTypes::StrainTensor &e, 
                          typename TensorTypes::BTensor &B, 
			  typename TensorTypes::DBTensor &DB,
			  typename TensorTypes::GradUTensor &gradU, 
			  typename TensorTypes::GradUDerivTensor &dgradUdqk) = 0;
  virtual void getE(typename TensorTypes::StrainTensor  &e, typename TensorTypes::GradUTensor &gradU) = 0;
  virtual bool isNonLinear() { return false; }
};


#endif
