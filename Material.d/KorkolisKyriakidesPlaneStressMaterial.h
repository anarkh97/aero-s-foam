// Sriramajayam

/*
 * KorkolisKyriakidesPlaneStressMaterial.h"
 * DG++
 *
 * Created by Ramsharan Rangarajan on 11/29/2010.
 *  
 * Copyright (c) 2006 Adrian Lew
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */ 

#ifndef KORKOLISKYRIAKIDESPLANESTRESSMATERIAL
#define KORKOLISKYRIAKIDESPLANESTRESSMATERIAL

#include <Material.d/ElastoPlasticPlaneStressMaterial.h>
#include <iostream>
#include <vector>
#include <cmath>


class KorkolisKyriakidesPlaneStressMaterial : public ElastoPlasticPlaneStressMaterial
{
 public:
  
  //! Constructor
  //! \param iLambda Lame constant for elastic response
  //! \param iMu Lame constant for elastic response
  //! \param iSigmaY Yield stress in 1D
  //! \param iK Isotropic hardening modulus
  //! \param iH Kinematic hardening modulus
  KorkolisKyriakidesPlaneStressMaterial(double iLambda, double iMu, 
					double iSigmaY, double iK = 0., 
					double iH = 0.);
  
  //! Destructor
  virtual ~KorkolisKyriakidesPlaneStressMaterial();

  //! Copy constructor
  //! \param MatObj Object to be copied
  KorkolisKyriakidesPlaneStressMaterial(const KorkolisKyriakidesPlaneStressMaterial &MatObj);

  //! Cloning
  virtual KorkolisKyriakidesPlaneStressMaterial * Clone() const;
  
  //! Compute the elastoplastic constitutive response
  //! Returns true if calculations went well and false otherwise
  //! \param Fnp1 Input. Deformation gradient at new state of material. Size 9x1
  //! \param CauchyStress Output. Has size 9x1
  //! \param Cep Output. Algorithmic elastoplastic tangent. If requested, has size 81x1
  //! \param UpdateFlag Input. Material state updated if true. Set to true by default.
  virtual bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1, 
							std::vector<double> * CauchyStress, 
							std::vector<double> * Cep = 0,
							const bool UpdateFlag = true);
  
  //! Returns the plastic strain in material (9x1 vector)
  virtual std::vector<double> GetMaterialPlasticStrain() const;
  
  //! Returns the equivalent plastic strain in material
  virtual double GetMaterialEquivalentPlasticStrain() const;
  
  //! Returns back stress in material (9x1 vector)
  virtual std::vector<double> GetMaterialBackStress() const;

  //! Returns the Isotropic hardening modulus
  virtual double GetIsotropicHardeningModulus() const;
  
  //! Returns the kinematic hardening modulus
  virtual double GetKinematicHardeningModulus() const;
  
  //! Returns the yield stress in 1D tension test
  virtual double GetYieldStressFromTensionTest() const;

  //! Returns the bulk modulus of material
  virtual double GetBulkModulus() const;
  
  //! Returns shear modulus of material
  virtual double GetShearModulus() const;
  
  //! Checks if the state of the material lies within the yield surface.
  //! \param CS Input. Cauchy stress 9x1 vector
  //! \param TOL Input. Tolerance to use for check
  //! The tolerance is non-dimensional. The check performed is 
  //! \f[\frac{f}{\sigma_Y}<TOL~\Rightarrow~\text{material state OK}. \f]
  virtual bool CheckMaterialState(const std::vector<double> &CS, const double TOL = 1.e-6) const;

 protected:
  
  // Compute the elastic constitutive response of material
  //! Returns true if calculations went well and false otherwise
  //! \param EPS Input. Elastic strain, 3x1
  //! \param CS Output. Computed cauchy stress, 3x1
  //! \param C Ouput. Elastic modulii 3x3
  virtual bool ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
						  std::vector<double> *CS, 
						  std::vector<double> *C=0) const;
  
  //! Evaluates the yield function
  //! \param Xi Input.  \f$\xi = \sigma-\sigma^b\f$. Size 3x1.
  //! \param eqP input. Equivalent plastic strain.
  virtual double EvaluateYieldFunction(const double * Xi, const double eqP) const;

  
 private:
  //! Evaluate derivative of yield function with respect to cauchy
  //! stress given \f$\xi\f$ and \f$e^p\f$.
  //! The derivative is computed numerically with a central difference formula.
  //! The perturbation equals \f$1.e-6 \times \sigma_Y\f$.
  //! \param Xi Input. \f$\xi\f$, size 3x1.
  //! \param eqP Input. \f$e^p\f$.
  std::vector<double> EvaluateDerivativeOfYieldFunction(const double * Xi, const double eqP) const;
 
  //! Computes the two linear transformations of stress involved in
  //! this model.
  //! \param Xi Input. \f$\xi\f$, size 3x1.
  //! \param L1Xi Output. First transformation of \f$\xi\f$, size 3x1.
  //! \param L2Xi Output. Second transformation of \f$\xi\f$, size 3x1.
  void ComputeLinearTransformationsOfStress(const double * Xi, 
					    double * L1Xi, double * L2Xi) const;
  
  //! Computes the two eigenvalues of a symmetric 2x2 tensr given in vector form.
  //! \param X Input. Vector form of symmetric tensor, size 3x1
  //! \param A Output. First eigenvalue
  //! \param B Output. Second eigenvalue
  void ComputeEigenvalues(const double * X, double &A, double &B) const;

  //! Computes the tensor norm of a deviatoric stress tensor
  //! given in vector form
  //! \param S Input. Symmetric deviator stress tensor, size 3x1.
  double DeviatoricStressNorm(const double * S) const;
  
  //! Computes the tensor norm of a deviatoric strain tensor
  //! given in vector form
  //! \param S Input. Symmetric deviator strain tensor, size 3x1.
  double DeviatoricStrainNorm(const double * S) const;
  
  
  //! Associative flow rule to compute \f$\xi\f$, given
  //! \f$\Delta \lambda\f$ and \f${\bf N}\f$.
  //! \param Xitrial Input. \f$\xi_{trial}\f$, size 3x1
  //! \param lambda Input. Value of \f$\Delta \lambda\f$.
  //! \param N Input. Value of \f${\bf N}\f$, size 3x1
  //! \param Xi output. Computed value of \f$\xi\f$.
 void ComputeXiGivenConsistencyParameterAndDirection(const double * Xitrial, const double lambda, 
						     const double * N, double * Xi) const;

 //! Associative flow rule to compute equivalent plastic strain 
 //! \f$e^p\f$ given \f$\delta \lambda\f$ and \f${\bf N}\f$.
 //! \param lambda Input. Value of \f$\Delta \lambda\f$.
 //! \param N Input. Value of \f${\bf N}\f$, size 3x1
 double ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(const double lambda, 
									    const double * N) const;
 
 
 //! Finds \f$\Delta\lambda\f$ such that \f$f=0\f$, given a direction \f${\bf N}\f$
 //! for the evolution of strain.
 //! Returns true if solution was found and false otherwise.
 //! \param Xitrial Input. \f$\xi_{trial}\f$, size 3x1
 //! \param N Input. Value of \f${\bf N}\f$, size 3x1
 //! \param lambda_L Input. Guess for \f$\Delta\lambda\f$ for which \f$f\leq 0\f$.
 //! \param lambda_R Input. Guess for \f$\Delta\lambda\f$ for which \f$f\geq 0\f$.
 //! \param lambda Output. Solved value of \f$\Delta\lambda\f$.
 bool ComputeConsistencyParameterGivenDirection(const double * Xitrial, const double * N, 
						double lambda_L, double lambda_R, double &lambda) const;

 
  //! Young's modulus
  double E;
  
  //! Poisson ratio
  double nu;
  
  //! Yield stress in 1D
  double SigmaY;
  
  //! Isotropic hardening modulus
  double K;
  
  //! Kinematic hardening modulus
  double H;
  
  //! Plastic strain 
  std::vector<double> EPSplastic;
  
  //! Back stress
  std::vector<double> BackStress;
  
  //! Equivalent plastic strain
  double equivEPSplastic;
};

#endif