#ifndef _STRAINEVALUATOR_H_
#define _STRAINEVALUATOR_H_
#include <Utils.d/NodeSpaceArray.h>

//Computes the geometrical part of the stiffness matrix 

class StrainEvaluator
{
  public:
    virtual Tensor *getTMInstance() = 0;
    virtual Tensor *getStressInstance() = 0;
    virtual Tensor *getStrainInstance() = 0;
    virtual Tensor *getBInstance(int numdofs) = 0;
    virtual Tensor *getDBInstance(int numdofs) = 0;
    virtual Tensor *getTempInstance(int numdofs) { return NULL; }
    virtual void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp) = 0;
    virtual void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp) = 0;
    virtual void getE(Tensor &e, Tensor &gradU) = 0;
};

class LinearStrain : public StrainEvaluator
{
  // to be used when the appropriate strain measure is the
  // infintesimal strain tensor e = 0.5*(gradU^t + gradU)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};

class GreenLagrangeStrain : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the
  // Green-Lagrange strain tensor E = 0.5*(gradU^t + gradU + gradU^T gradU) which is a symmetric rank 2 tensor
  // Constitutive models based on this strain measure should return 
  // the second Piola-Kirchoff stress tensor S (rank 2, symmetric) and
  // the second (material) elasticity tensor M (rank 4, major and minor symmetries)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getTempInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};

class LogarithmicStrain : public StrainEvaluator
{
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};

class PrincipalStretches : public StrainEvaluator
{
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};

class LogarithmicPrincipalStretches : public StrainEvaluator
{
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};

class DeformationGradient : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the
  // deformation gradient tensor F = I + gradU, which is a nonsymmetric rank 2 tensor
  // Constitutive models based on this strain measure should return 
  // the first Piola-Kirchoff stress tensor P (rank 2, nonsymmetric) XXXX or should it be the nominal stress tensor note: P^T = N !!!!!
  // and the first elasticity tensor A (rank 4, major symmetries only)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk, Tensor *temp);
    void getE(Tensor &e, Tensor &gradU);
};


template<class TensorTypes>
class GenStrainEvaluator
{
  public:
    virtual void getEBandDB(typename TensorTypes::StrainTensor &e, 
                            typename TensorTypes::BTensor &B, 
                            typename TensorTypes::DBTensor &DB,
                            typename TensorTypes::GradUTensor &gradU, 
                            typename TensorTypes::GradUDerivTensor &dgradUdqk) = 0;
    virtual void getEandB(typename TensorTypes::StrainTensor &e,
                          typename TensorTypes::BTensor &B,
                          typename TensorTypes::GradUTensor &gradU,
                          typename TensorTypes::GradUDerivTensor &dgradUdqk) = 0;
    virtual void getE(typename TensorTypes::StrainTensor &e, typename TensorTypes::GradUTensor &gradU) = 0;
    virtual bool isNonLinear() { return false; }
};

#endif
