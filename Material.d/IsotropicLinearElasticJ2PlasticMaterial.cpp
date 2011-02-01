// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticMaterial.cpp
 * DG++
 *
 * Created by Ramsharan Rangarajan on 09/21/2010.
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


#include "IsotropicLinearElasticJ2PlasticMaterial.h"


// Constructor
IsotropicLinearElasticJ2PlasticMaterial::
IsotropicLinearElasticJ2PlasticMaterial(double iLambda, double iMu, 
					double iSigmaY, double iK, 
					double iH) 
  :Mu(iMu), SigmaY(iSigmaY), K(iK), H(iH)
{
  // Compute bulk modulus
  Kappa = iLambda + (2./3.)*iMu;
  
  // Create material for elastic response
  ILE = new IsotropicLinearElastic(iLambda, iMu);
  
  // Initialize plastic strain to zero.
  EPSplastic.clear();
  for(int i=0; i<9; i++)
    EPSplastic.push_back( 0. );
  
  // Initialize equivalent plastic strain to zero
  equivEPSplastic = 0.;
  
  // Initialize back stress to zero
  BackStress.clear();
  for(int i=0; i<9; i++)
    BackStress.push_back( 0. );
}


// Destructor
IsotropicLinearElasticJ2PlasticMaterial::
~IsotropicLinearElasticJ2PlasticMaterial()
{ delete ILE; }


// Copy constructor
IsotropicLinearElasticJ2PlasticMaterial::
IsotropicLinearElasticJ2PlasticMaterial(const IsotropicLinearElasticJ2PlasticMaterial &Mat)
  : Kappa(Mat.Kappa), Mu(Mat.Mu), SigmaY(Mat.SigmaY), K(Mat.K), H(Mat.H)
{ 
  ILE = new IsotropicLinearElastic(*Mat.ILE); 
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<9; i++)
    {
      EPSplastic.push_back( Mat.EPSplastic[i] );
      BackStress.push_back( Mat.BackStress[i] );
    }
  equivEPSplastic = Mat.equivEPSplastic;
}


// Cloning
IsotropicLinearElasticJ2PlasticMaterial *
IsotropicLinearElasticJ2PlasticMaterial::Clone() const
{ return new IsotropicLinearElasticJ2PlasticMaterial(*this); }

// Return plastic strain in material
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::GetMaterialPlasticStrain() const
{ return EPSplastic; }

// Resturn the equivalent plastic strain in material
double IsotropicLinearElasticJ2PlasticMaterial::GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }

// Return back stress in material
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::GetMaterialBackStress() const
{ return BackStress; }


// Return isotropic hardening modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetIsotropicHardeningModulus() const
{ return K; }

// Return kinematic hardening modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetKinematicHardeningModulus() const
{ return H; }

// Return flow stress from uniaxial tension test
double IsotropicLinearElasticJ2PlasticMaterial::GetYieldStressFromTensionTest() const
{ return SigmaY; }


// Return bulk modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetBulkModulus() const
{ return Kappa; }


// Return shear modulus of material
double IsotropicLinearElasticJ2PlasticMaterial::GetShearModulus() const
{ return Mu; }


// Compute linear elastic response
bool IsotropicLinearElasticJ2PlasticMaterial::
ComputeElasticConstitutiveResponse(const std::vector<double> &EPS, 
				   std::vector<double> *CS, 
				   std::vector<double> *C) const
{ 
  // Resize of necessary
  if( int(CS->size())<9 )
    CS->resize( 9 );

  if( C )
    if( int(C->size())<81 )
      C->resize( 81 );
  
  // Convert strain to deformation gradient
  std::vector<double> F = EPS;
  F[0] += 1.;
  F[4] += 1.;
  F[8] += 1.;

  return ILE->GetConstitutiveResponse(&F, CS, C); 
}


// Compute deviatoric part of tensor
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::Deviatoric(const double * T) const
{
  // 1/3 of trace of T
  double TRby3 = (T[0] + T[4] + T[8])/3.;

  // Identity
  double I[] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
 
  // Compute deviatoric part of T as T-(1/3)trace(T)*I.
  std::vector<double> devT(9);
  for(int i=0; i<9; i++)
    devT[i] = T[i]-TRby3*I[i];
  
  return devT;
}


// Compute norm of tensor
double IsotropicLinearElasticJ2PlasticMaterial::Norm(const double *T) const
{ 
  double norm2 = 0.;
  for(int i=0; i<9; i++)
    norm2 += pow(T[i],2.);
  return sqrt(norm2);
}


// Evaluate yield function
double IsotropicLinearElasticJ2PlasticMaterial::
EvaluateYieldFunction(const double * CauchyStress,
		      const double * SigmaB, 
		      const double eqP) const
{
  // f = norm(xi) - sqrt(2/3) * (SigmaY + K*eqP),
  // where 
  // xi = dev(Cauchy stress)-SigmaB
  // eqP = equivalent plastic strain 
  
  // Deviatoric part of CauchyStress
  std::vector<double> S = Deviatoric(&CauchyStress[0]);
  
  // xi = S-SigmaB
  double xi[9];
  for(int i=0; i<9; i++)
    xi[i] = S[i]-SigmaB[i];

  // Norm of xi
  double normXI = Norm(xi);
  
  return normXI - sqrt(2./3.)*(SigmaY+K*eqP);
}
  
  
// Compute elasto-plastic constitutive response
bool IsotropicLinearElasticJ2PlasticMaterial::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1, 
					 std::vector<double> *CauchyStress, 
					 std::vector<double> *Cep, 
					 const bool UpdateFlag)
{
  // Resize output cauchy stress
  if( int(CauchyStress->size())<9 )
    CauchyStress->resize(9);
  
  // Resize output tangents if requested
  if( Cep )
    if( int(Cep->size())<81 )
      Cep->resize( 81 );
  
  // Elastic modulii
  std::vector<double> * Ce = 0;
  if( Cep )
    Ce = new std::vector<double>(81);

  // Identity tensor
  double I[9] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
  
  // Compute symmetric infinitesimal strain tensor
  std::vector<double> EPS(9);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      EPS[3*i+j] = 0.5*(Fnp1[3*i+j]+Fnp1[3*j+i])-I[3*i+j];
  
  // EVALUATE TRIAL STATE BY FREEZING PLASTIC FLOW
  
  // Trial elastic strain
  std::vector<double> epsEtrial(9);
  for(int i=0; i<9; i++)
    epsEtrial[i] = EPS[i]-EPSplastic[i];
  
  // Trial cauchy stress
  std::vector<double> CStrial(9);
  if( !ComputeElasticConstitutiveResponse(epsEtrial, &CStrial, Ce) )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
	       <<"Could not compute elastic response of material. \n";
      return false;
    }
  
  // Evaluate yield function at trial state
  double Ftrial = EvaluateYieldFunction(&CStrial[0], &BackStress[0], equivEPSplastic);
  
  
  // COMPUTE NEW STATE OF MATERIAL
  
  if( Ftrial<= 0.)
    {
      // This step is purely elastic.
      // Final state of the material is the trial state.
      
      for(int i=0; i<9; i++)
	(*CauchyStress)[i] = CStrial[i];
      
      if( Cep )
	for(int i=0; i<81; i++)
	  (*Cep)[i] = (*Ce)[i];
      
      // No need to update plastic internal variables
    }
  
  else
    {
      // This step is elasto-plastic.
      // The trial state is not the actual state of the material
      
      // Compute consistency parameter
      double dLambda = Ftrial/(2.*Mu + (2./3.)*(H+K));
      
      // Strial
      std::vector<double> Strial = Deviatoric(&CStrial[0]);
      
      // XItrial = Strial-BackStress
      std::vector<double> XItrial(9);
      for(int i=0; i<9; i++)
	XItrial[i] = Strial[i]-BackStress[i];
            
      double normXItrial = Norm(&XItrial[0]);
      if( normXItrial<1.e-8 )
	{
	  std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
		   <<"Norm of trial stress close to zero. \n";
	  return false;
	}
      
      // Unit tensor XItrial/norm(XItrial) = XI_(n+1)/norm(X_(n+1)).
      double N[9];
      for(int i=0; i<9; i++)
	N[i] = XItrial[i]/normXItrial;
      
      // New elastic strain
      std::vector<double> NewEPSelastic(9);
      
      for(int i=0; i<9; i++)
	NewEPSelastic[i] = EPS[i] - EPSplastic[i] - dLambda*N[i];
      
      // New Cauchy Stress
      if( !ComputeElasticConstitutiveResponse(NewEPSelastic, CauchyStress) )
	{
	  std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
		   <<"Could not compute elastic response.\n";
	  return false;
	}
      
      if( Cep )
	{
	  // Evaluate the consistency eplasto-plastic modulii
	  double theta    = 1. -2.*Mu*dLambda/normXItrial;
	  double thetabar = 1./(1.+(K+H)/(3.*Mu)) + theta -1.;
	  
	  for(int i=0; i<3; i++)
	    for(int j=0; j<3; j++)
	      for(int k=0; k<3; k++)
		for(int L=0; L<3; L++)
		  (*Cep)[27*i+9*j+3*k+L] = 
		    Kappa*I[3*i+j]*I[3*k+L] + 
		    2.*Mu*theta*( 0.5*(I[3*i+k]*I[3*j+L] + I[3*i+L]*I[3*j+k]) - (1./3.)*I[3*i+j]*I[3*k+L] ) -
		    2.*Mu*thetabar*N[3*i+j]*N[3*k+L];
	}
	  
      // Update plastic strain and backstress if requested
      if( UpdateFlag==true )
	{
	  equivEPSplastic += sqrt(2./3.)*dLambda;
	  for(int i=0; i<9; i++)
	    {
	      EPSplastic[i] += dLambda*N[i];
	      BackStress[i] += (2./3.)*H*dLambda*N[i];
	    }
	}
    }
  
  if( Ce )
    delete Ce;
  
  return true;
}


bool IsotropicLinearElasticJ2PlasticMaterial::
CheckMaterialState(const std::vector<double> &CS, const double TOL) const
{
  double f = EvaluateYieldFunction(&CS[0], &BackStress[0], equivEPSplastic);
  if( f/SigmaY < TOL )
    return true;
  else
    return false;
}
  

