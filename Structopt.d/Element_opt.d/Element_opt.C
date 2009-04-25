#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Element_opt.h>

void
Element_opt::getGradVonMises(Vector &dstress, Vector &weight, 
			     CoordSet &cs, CoordSet &dcs, 
			     Vector &elDisp, Vector &elGrad,
			     int strInd, int surface, 
			     double*, double*)
{
  dstress.zero();
  weight.zero();
}

void
Element_opt::getGradVonMisesInt(CoordSet& cs, CoordSet& dcs, Vector& d, Vector& dd,
				double& sigbar,double& fac,int areaFlag,
				double& vmint,double& dvmint,double& vol,double& dvol,
				double* , double*)
{
  vmint  = 0.0;
  vol    = 0.0;
  dvmint = 0.0;
  dvol   = 0.0;
}

void
Element_opt::computeDfDsPressure(CoordSet&, CoordSet&, Vector& eleDfDs, 
				 GeomState *, int cflg)
{
  eleDfDs.zero();
}

void
Element_opt::computeDfDsThermal(CoordSet&, CoordSet&, Vector &, Vector &,
				Vector& eleDfDs,int glflag, GeomState *,int)
{
  eleDfDs.zero();
}

//------------------------------------------------------------------------------
void Element_opt::getGradElectricField(CoordSet &cs,CoordSet &dcs,
				       Vector &gradEfieldx,Vector &gradEfieldy,
				       Vector &gradEfieldz,Vector &weight,
				       Vector &d_weight,Vector &volt,
				       Vector &gradVolt)
{
  
  weight.zero();
  d_weight.zero();
  
  gradEfieldx.zero();
  gradEfieldy.zero();
  gradEfieldz.zero();
  
}

//------------------------------------------------------------------------------

void Element_opt::getGradElectricForce(CoordSet &cs,CoordSet &dcs,Vector &efieldx,
				       Vector &efieldy,Vector &efieldz,
				       Vector &gradEfieldx,Vector &gradEfieldy,
				       Vector &gradEfieldz,Vector &gradEforcex,
				       Vector &gradEforcey,Vector &gradEforcez)
{
  gradEforcex.zero();
  gradEforcey.zero();
  gradEforcez.zero();
}

//------------------------------------------------------------------------------

void Element_opt::computeGradSourceFlux(Vector &d_elForce,CoordSet&,Vector&,
					Vector&, Vector&)
{
  d_elForce.zero();
}


void Element_opt::gradstiffness(CoordSet &, CoordSet &, FullSquareMatrix& result, int) 
{ 
  fprintf(stderr,"WARNING: using zero stiffness matrix\n"); 
  result.setSize(numDofs());
  result.zero();  
} 

void Element_opt::gradMassMatrix(CoordSet &, CoordSet &, FullSquareMatrix& result, double mratio)
{ 
  fprintf(stderr,"WARNING: using zero mass matrix\n");
  result.setSize(numDofs());
  result.zero();
}


void Element_opt::gradDampMatrix(CoordSet &, CoordSet &, FullSquareMatrix& result,  double freq) 
{ 
  //fprintf(stderr,"WARNING: using zero damping matrix\n");
  result.setSize(numDofs());
  result.zero();
}


void Element_opt::getGradVonMises(ComplexVector &dstress, Vector &weight, 
				  CoordSet &cs, CoordSet &dcs, 
				  ComplexVector &elDisp, ComplexVector &elGrad,
				  int strInd, int surface, 
				  double* ndTemp, double* dndTemp)
{
  Vector relDisp(elDisp.size());
  Vector relGrad(elGrad.size());
  Vector rdStress(dstress.size());  
  for(int i=0; i<relDisp.size(); ++i)
    { 
      relDisp[i] = elDisp[i].real(); 
      relGrad[i] = elGrad[i].real(); 
    }
  getGradVonMises(rdStress, weight, cs, dcs, relDisp, relGrad, strInd, surface, ndTemp, dndTemp);
  for(int i=0; i<relDisp.size(); ++i)
    { dstress[i] = DComplex(rdStress[i], 0.0); }
  for(int i=0; i<relDisp.size(); ++i)
    { 
      relDisp[i] = elDisp[i].imag(); 
      relGrad[i] = elGrad[i].imag(); 
    }
  getGradVonMises(rdStress, weight, cs, dcs, relDisp, relGrad, strInd, surface, ndTemp, dndTemp);
  for(int i=0; i<relDisp.size(); ++i)
    { dstress[i] += DComplex(0.0, rdStress[i]); }
  return;
}

//------------------------------------------------------------------------------
void Element_opt::rotateConstitutiveMatrix33(double Cin[3][3], double *T33, double Cout[3][3])
{
  double l1 = T33[0]; double m1 = T33[1]; double n1 = T33[2];
  double l2 = T33[3]; double m2 = T33[4]; double n2 = T33[5];
  double l3 = T33[6]; double m3 = T33[7]; double n3 = T33[8];

  double TNN[3][3] = {{l1*l1,m1*m1,n1*n1},{l2*l2,m2*m2,n2*n2},{l3*l3,m3*m3,n3*n3}};
  double CC[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};  
  // Cout = T'(CT)  
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<3; ++k) 
	    {
	      CC[i][j] += Cin[i][k]*TNN[k][j];
	    }
	}
    }
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  Cout[i][j] = 0.0;
	  for(int k=0; k<3; ++k) 
	    {
	      Cout[i][j] += TNN[k][i]*CC[k][j];
	    }
	}
    }
  return;
}

//------------------------------------------------------------------------------
void Element_opt::rotateConstitutiveMatrix63(double Cin[6][3], double *T33, double Cout[6][3])
{
  double l1 = T33[0]; double m1 = T33[1]; double n1 = T33[2];
  double l2 = T33[3]; double m2 = T33[4]; double n2 = T33[5];
  double l3 = T33[6]; double m3 = T33[7]; double n3 = T33[8];
  
  double TNN[3][3] = {{l1*l1,m1*m1,n1*n1},{l2*l2,m2*m2,n2*n2},{l3*l3,m3*m3,n3*n3}};
  double TNS[3][3] = {{l1*m1,l1*n1,m1*n1},{l2*m2,l2*n2,m2*n2},{l3*m3,l3*n3,m3*n3}};
  double TSN[3][3] = {{l1*l2,m1*m2,n1*n2},{l1*l3,m1*m3,n1*n3},{l2*l3,m2*m3,n2*n3}};
  double TSS[3][3] = {{l1*m2+l2*m1,l1*n2+l2*n1,m1*n2+m2*n1},
                      {l1*m3+l3*m1,l1*n3+l3*n1,m1*n3+m3*n1},
                      {l2*m3+l3*m2,l2*n3+l3*n2,m2*n3+m3*n2}};
  double CC[6][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};  

  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{ 
	  TNS[i][j] *= M_SQRT2;  
	  TSN[i][j] *= M_SQRT2; 
	}
    }

  double T66[6][6] = {{TNN[0][0],TNN[0][1],TNN[0][2],TNS[0][0],TNS[0][1],TNS[0][2]},
                      {TNN[1][0],TNN[1][1],TNN[1][2],TNS[1][0],TNS[1][1],TNS[1][2]},  
                      {TNN[2][0],TNN[2][1],TNN[2][2],TNS[2][0],TNS[2][1],TNS[2][2]},  
                      {TSN[0][0],TSN[0][1],TSN[0][2],TSS[0][0],TSS[0][1],TSS[0][2]},
                      {TSN[1][0],TSN[1][1],TSN[1][2],TSS[1][0],TSS[1][1],TSS[1][2]},  
                      {TSN[2][0],TSN[2][1],TSN[2][2],TSS[2][0],TSS[2][1],TSS[2][2]}};
  
  for(int j=3; j<6; ++j)
    {
      for(int i=0; i<3; ++i) 
	{ Cin[j][i] *= M_SQRT2; }
    }
  for(int i=0; i<6; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<3; ++k)
	    {
	      CC[i][j] += Cin[i][k]*TNN[k][j];
	    }
	}
    }
  for(int i=0; i<6; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  Cout[i][j] = 0.0;
	  for(int k=0; k<6; ++k)
	    {
	      Cout[i][j] += T66[i][k]*CC[k][j];
	    }
	}
    }
  for(int j=3; j<6; ++j)
    {
      for(int i=0; i<3; ++i) 
	{ Cout[j][i] *= M_SQRT1_2; }
    }
  return;
}

/*
//------------------------------------------------------------------------------
void Element_opt::gradstiffness(CoordSet & cs, CoordSet & dcs, FullSquareMatrixC& result, int cmd) 
{
  double* dstiff_d = static_cast<double*>(dbg_alloca(sizeof(double)*numDofs()*numDofs()));
  FullSquareMatrix result_d(1, dstiff_d);
  result_d.setSize(numDofs());
  gradstiffness(cs, dcs, result_d, cmd);
  result = result_d;
} 

void Element_opt::gradMassMatrix(CoordSet & cs, CoordSet & dcs, FullSquareMatrixC& result, double mratio)
{ 
  double* dmass_d = static_cast<double*>(dbg_alloca(sizeof(double)*numDofs()*numDofs()));
  FullSquareMatrix result_d(1, dmass_d);
  result_d.setSize(numDofs());
  gradMassMatrix(cs, dcs, result_d);
  result = result_d;
}

void Element_opt::gradDampMatrix(CoordSet & cs, CoordSet & dcs, FullSquareMatrixC& result, double freq) 
{ 
  double* ddamp_d = static_cast<double*>(dbg_alloca(sizeof(double)*numDofs()*numDofs()));
  FullSquareMatrix result_d(1, ddamp_d);
  result_d.setSize(numDofs());
  gradDampMatrix(cs, dcs, result_d, freq);
  result = result_d;
}
*/

#endif
