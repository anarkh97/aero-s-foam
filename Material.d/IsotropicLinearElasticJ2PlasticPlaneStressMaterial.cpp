// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticPlaneStressMaterial.cpp
 * DG++
 *
 * Created by Ramsharan Rangarajan on 11/08/2010.
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

#include "IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h"


// Constructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
IsotropicLinearElasticJ2PlasticPlaneStressMaterial(double iLambda, double iMu, 
						   double iSigmaY, double iK,
						   double iH)
{
  // Youngs modulus
  E = iMu*(3.*iLambda + 2.*iMu)/(iLambda + iMu);
  
  // Poisson ratio
  nu = 0.5*iLambda/(iLambda+iMu);

  // Constants related to plasticity
  SigmaY = iSigmaY;
  K = iK;
  H = iH;
  
  // Zero initial plastic strain and backstress
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++)
    {
      EPSplastic.push_back( 0. );
      BackStress.push_back( 0. );
    }
  equivEPSplastic = 0.;
}


// Destructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::~IsotropicLinearElasticJ2PlasticPlaneStressMaterial()
{}


// Copy constructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
IsotropicLinearElasticJ2PlasticPlaneStressMaterial(const IsotropicLinearElasticJ2PlasticPlaneStressMaterial &Mat)
  :E(Mat.E), nu(Mat.nu), SigmaY(Mat.SigmaY), K(Mat.K), H(Mat.H)
{
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++)
    {
      EPSplastic.push_back( Mat.EPSplastic[i] );
      BackStress.push_back( Mat.BackStress[i] );
    }
  equivEPSplastic = Mat.equivEPSplastic;
}



// Cloning
IsotropicLinearElasticJ2PlasticPlaneStressMaterial * 
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::Clone() const
{ return new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(*this); }


// Return plastic strain
std::vector<double> IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialPlasticStrain() const
{
  std::vector<double> EP(9,0.);
  EP[0] = EPSplastic[0];
  EP[4] = EPSplastic[1];
  EP[8] = -(EP[0]+EP[4]);
  EP[1] = EP[3] = 0.5*EPSplastic[2];
  EP[2] = EP[5] = EP[6] = EP[7] = 0.;
  return EP;
}


// Return equivalent plastic strain
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }


// Return back stress
std::vector<double> IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialBackStress() const
{
  std::vector<double> BS(9,0.);
  BS[0] = BackStress[0];
  BS[4] = BackStress[1];
  BS[8] = -(BS[0]+BS[4]);
  BS[1] = BS[3] = BackStress[2];
  BS[2] = BS[5] = BS[6] = BS[7] = 0.;
  return BS;
}


// Return isotropic hardening modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetIsotropicHardeningModulus() const
{ return K; }

// Return kinematic hardening modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetKinematicHardeningModulus() const
{ return H; }

// Return yield stress for 1D tension test
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetYieldStressFromTensionTest() const
{ return SigmaY; }

// Return Bulk modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetBulkModulus() const
{ return E/(3.*(1.-2.*nu)); }

// Return shear modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetShearModulus() const
{ return 0.5*E/(1.+nu); }


// Compute the elastic constitutive response
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
ComputeElasticConstitutiveResponse(const std::vector<double> &EPS, 
				   std::vector<double> *CS, 
				   std::vector<double> *C) const
{
  if( int(CS->size())<3 )
    CS->resize(3);
  
  double v = E/(1.-nu*nu);
  double Cmat[3][3] = {{v*1., v*nu, 0.},
		       {v*nu, v, 0.},
		       {0., 0., v*(1.-nu)/2.}};
  
  for(int i=0; i<3; i++)
    {
      (*CS)[i] = 0.;
      for(int j=0; j<3; j++)
	(*CS)[i] += Cmat[i][j]*EPS[j];
    }
  
  if( C!=0 )
    {
      if( int(C->size())<9 )
	C->resize( 9 );
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  (*C)[3*i+j] = Cmat[i][j];
    }
 
  return true;
}


// Evaluate norm of deviatoric part of \f$\sigma-\sigma^b\f$
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeJ2(const double *Xi) const
{
  double P[3][3] = {{2./3.,  -1./3., 0.},
		    {-1./3., 2./3.,  0.},
		    {0., 0., 2.}};
  
  // Norm of deviatoric part of Xi = Transpose(Xi) * P * Xi.
  double normXi2 = 0.;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      normXi2 += Xi[i]*P[i][j]*Xi[j];
  return sqrt(normXi2);
}


// Evaluate yield function
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
EvaluateYieldFunction(const double * Xi,  const double eqP) const
{ 
  // Norm of Xi
  double J2 = ComputeJ2(Xi);
  
  // Evaluate radius of yield surface
  double YSrad = sqrt(2./3.) * (SigmaY + K*eqP);
  
  return J2 - YSrad;
}




// Check if state of material lies within yield surface
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
CheckMaterialState(const std::vector<double> &CS, 
		   const double TOL) const
{
  // Compute Xi = Cauchy stress - Back stress (in vector form)
  double Xi[3] = { CS[0]-BackStress[0],
		   CS[4]-BackStress[1],
		   CS[1]-BackStress[2] };
  
  // Evaluate yield function
  double Fval = EvaluateYieldFunction(Xi, equivEPSplastic);
  
  // Check for tolerance after non-dimensionalizing Fval.
  if( Fval/SigmaY < TOL )
    return true;
  else
    return false;
}


// Compute elastoplastic response of material
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
					 std::vector<double> *CauchyStress,
					 std::vector<double> *Cep,
					 const bool UpdateFlag) 
{
  // Notify that elastoplastic tangents are not computed.
  if( Cep!=0 )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
	       <<"Elastoplastic tangents not implemented.\n";
      return false;
    }
  
  
  // Vector form of symmetric strain
  double EPS[3] = { Fnp1[0]-1., Fnp1[4]-1., Fnp1[1]+Fnp1[3] };
  
  // Compute trial elastic state by freezing plastic flow
  std::vector<double> EPSelas(3,0.);
  for(int i=0; i<3; i++)
    EPSelas[i] = EPS[i]-EPSplastic[i];
  
  // Compute trial cauchy stress
  std::vector<double> CStrial(3,0.);
  
  if( !ComputeElasticConstitutiveResponse(EPSelas, &CStrial) )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastioPlasticConstitutiveResponse()- "
	       <<"Could not compute elastic response at trial state.\n";
      return false;
    }
  
  // Trial Xi
  double Xitrial[] = { CStrial[0]-BackStress[0],
		       CStrial[1]-BackStress[1],
		       CStrial[2]-BackStress[2] };
  
  // Evaluate yield function at trial state
  double Ftrial = EvaluateYieldFunction(Xitrial, equivEPSplastic);
  
  // Use some tolerance for checking yield function value
  double TOL = SigmaY*1.e-6;
  
  // Resize outputs if required
  if( int(CauchyStress->size())<9 )
    CauchyStress->resize(9);
  
  if( Ftrial<TOL )
    {
      // This step is purely elastic.
      // No need to update plastic variables.
      
      // Convert cauchy stress to tensor form
      (*CauchyStress)[0] = CStrial[0];
      (*CauchyStress)[4] = CStrial[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CStrial[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;
      
      return true;
    }
  
  else
    {
      // This step is elasto-plastic.
      // Need to compute consistency parameter, called lambda here.
      // Use bisection method to compute lambda.
      
      // Max. number of iteration = 1000.
      int nItMax = 1000;
      // Tolerance for checking convergence if f<TOL
      bool CONVERGED = false;
      
      // lambda = 0 is the trial state and Ftrial is surely > 0.
      double lambda_R = 0.;
      double F_R = Ftrial;
      
      // The other extreme is when there is plenty of plastic strain.
      double epMAX = equivEPSplastic + sqrt(2./3.)*sqrt( pow(EPS[0]-EPSplastic[0],2.) + 
							 pow(EPS[1]-EPSplastic[1],2.) + 
							 0.5*pow(EPS[2]-EPSplastic[2],2.) );  // due the factor of 2 in the strains.
      
      // Step size for lambda to probe where F < 0.
      double lambdaStep = 1.5*(epMAX-equivEPSplastic)/(SigmaY+K*epMAX);  
      
      // For computing Xi as a function of lambda.
      double Xi[] = {0.,0.,0.};
      
      // Proble for lambda_L that corresponds to F < 0.
      double lambda_L = 0.;
      double F_L = 0.;
      bool NegFfound = false;
      int counter = 0;
      while( NegFfound==false )
	{
	  // Increment counter
	  counter++;
	  // Increment lamba by one step
	  lambda_L += lambdaStep;
	  
	  // Evaluate corresponding value of Xi
	  if( !ComputeXiGivenConsistencyParameter(Xitrial, lambda_L, Xi) )
	    {
	      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElasticConstitutiveResponse()- "
		       <<"Could not evaluate stress state for some guess of consistency parameter.\n";
	      return false;
	    }
	  
	  // Evaluate yield function
	  F_L = EvaluateYieldFunction(Xi, equivEPSplastic+sqrt(2./3.)*lambda_L*ComputeJ2(Xi));
	  
	  if( F_L<=0. )
	    // lambda corresponding to negative value of F has been found.
	    NegFfound = true;
	  else
	    {
	      // F is still positive.
	      // This means F = 0 lies to the left of lambda_L.
	      // So shift the right side of the bracket to this point.
	      lambda_R = lambda_L;
	      F_R      = F_L;
	    }
	  
	  // Sound a warning if lambda is orders of magnitude away from lambdaStep.
	  if( counter%100 == 0)
	    std::cerr<<"\n WARNING: "
		     <<"IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastioPlasticConstitutiveResponse()- "
		     <<" Consistency parameter exceeds "<<counter<<" times initial estimate.\n";
	  
	  // Upper limit on how far to probe- avoid an infinite loop.
	  if( counter > 1000 )
	    return false;
	}
      
      // Check that F has been bracketed.
      if( F_L>0. || F_R<0. )
	{
	  std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastioPlasticConstitutiveResponse()- "
		   <<"Could not bracket yield function to solve for consistency parameter.\n";
	  return false;
	}
      
      // Required solution lies between lambda_L and lambda_R
      double lambda = 0.;
      
      // Iterate
      for(int it=0; it<nItMax && !CONVERGED; it++)
	{
	  lambda = (lambda_L + lambda_R)/2.;
	  if( !ComputeXiGivenConsistencyParameter(Xitrial, lambda, Xi) )
	    {
	      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElasticConstitutiveResponse()- "
		       <<"Could not evaluate stress state for some guess of consistency parameter.\n";
	      return false;
	    }
	  double F = EvaluateYieldFunction(Xi, equivEPSplastic+sqrt(2./3.)*lambda*ComputeJ2(Xi));
	  
	  if( fabs(F) < TOL )
	    CONVERGED = true;
	  else
	    {
	      if( F<0. )
		{
		  lambda_L = lambda;
		  F_L      = F;
		}
	      else
		{
		  lambda_R = lambda;
		  F_R      = F;
		}
	    }
	}
      
      if( !CONVERGED )
	{
	  std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
		   <<"Iterations to compute consistency parameter did not converge.\n";
	  return false;
	}
      
      // Compute P*Xi
      double PXi[] = { (2./3.)*Xi[0] - (1./3.)*Xi[1], 
		       (2./3.)*Xi[1] - (1./3.)*Xi[0], 
		       2.*Xi[2] };
      
      // Compute new elastic strain
      std::vector<double> EPSelastic(3);
      for(int i=0; i<3; i++)
	EPSelastic[i] = EPS[i] - EPSplastic[i] - lambda*PXi[i];
      
      // Compute new cauchy stress
      std::vector<double> CS(3);
      if( !ComputeElasticConstitutiveResponse(EPSelastic, &CS, 0) )
	{
	  std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
		   <<"Could not compute elastic response.\n";
	  return false;
	}
      
      // Convert vector to tensor form
      (*CauchyStress)[0] = CS[0];
      (*CauchyStress)[4] = CS[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CS[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;
      
      // If requested, update state of material
      if( UpdateFlag==true )
	{
	  for(int i=0; i<3; i++)
	    {
	      // Update plastic strain
	      EPSplastic[i] += lambda*PXi[i];	
	      // Update back stress
	      BackStress[i] += (2./3.)*H*lambda*Xi[i];
	    }
	  // Update equivalent plastic strain
	  equivEPSplastic += sqrt(2./3.)*lambda*ComputeJ2(Xi);
	}

      return true;
    }
}



// Given a value of consistency parameter and trial Xi, compute new Xi
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
ComputeXiGivenConsistencyParameter(const double * Xitrial, const double lambda, double * Xi) const
{
  double v = lambda*E/(3.*(1.-nu*nu));
 
  double a = 1.+(2./3.)*lambda*H + v*(2.-nu);
  double b = v*(2.*nu-1.);
  double c = 1.+(2./3.)*lambda*H + v*3.*(1.-nu);
  
  // M*Xi = Xitrial,
  // where M = (1+(2/3)*lambda*H) + lambda*C*P
  // M = [a b 0]
  //     [b a 0]
  //     [0 0 c]
  
  double a2mb2 = pow(a,2.)-pow(b,2.);
  
  // Check that M is invertible.
  if( fabs(c*a2mb2) < 1.e-6 )
    return false;
  else
    {
      Xi[0] = (a*Xitrial[0]-b*Xitrial[1])/a2mb2;
      Xi[1] = (a*Xitrial[1]-b*Xitrial[0])/a2mb2;
      Xi[2] = Xitrial[2]/c;
      return true;
    }
}