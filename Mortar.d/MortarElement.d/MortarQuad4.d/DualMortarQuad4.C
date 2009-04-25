// ----------------------------------------------------------------
// HB - 08/13/03
// ----------------------------------------------------------------
// Standard C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// FEM headers
#include <Element.d/Element.h>
#include <Math.d/matrix.h>

#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarQuad4.d/DualMortarQuad4.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarQuad4::DualMortarQuad4() 
{
  Initialize();
} 

DualMortarQuad4::DualMortarQuad4(FaceElement* FaceElem)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
}

DualMortarQuad4::DualMortarQuad4(double _area, FaceElement* FaceElem)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}

DualMortarQuad4::DualMortarQuad4(FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
  //cerr << " In DualMortarQuad4::DualMortarQuad4: " << endl;
  //Alpha->print("Alpha=");
}

DualMortarQuad4::DualMortarQuad4(double _area, FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarQuad4::~DualMortarQuad4()
{
  if(Alpha) { delete Alpha; Alpha = 0; }
}

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS
// -----------------------------------------------------------------------------------------------------
void
DualMortarQuad4::Initialize()
{
  Alpha = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
DualMortarQuad4::nNodes() { return 4; }

int
DualMortarQuad4::nMortarShapeFct() { return 4; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
DualMortarQuad4::ComputeDualCoeffs(CoordSet &cs)
{
  if(Alpha==0) Alpha = new FullM(4);
  else{ delete Alpha; Alpha = new FullM(4); }
  Alpha->zero();

  // 1) compute scalar mass matrix (M) of supporting face element
  FaceElement* FaceElem = GetPtrMasterFace();
  FullM M = FaceElem->ScalarMass(cs, 1.0, 2);
  // 2) compute std shape function integrals (b) over supporting face element
  double ShapeIntg[4];
  FaceElem->IntegrateShapeFcts(ShapeIntg, cs, 1.0, 2);
  // 3) solve M.alpha=b to get the dual shape fcts coeffs alpha
  M.factor();
  for(int j=0; j<4; j++){
    (*Alpha)[j][j] = ShapeIntg[j];
    M.reSolve((*Alpha)[j]); 
  }
}

// !!! ASSUME THAT DUAL MORTAR COEFFS Alpha HAVE ALREADY BEEN COMPUTED !!!
void
DualMortarQuad4::GetDualMortarShapeFct(double* Shape, double* m)
{
   if(Shape==0) Shape = new double[4];
   double StdShape[4];
   GetPtrMasterFace()->GetShapeFctVal(StdShape, m);
   for(int i=0; i<4; ++i){
     Shape[i] = 0;
     for(int j=0; j<4; j++)
       Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
   }
}

void
DualMortarQuad4::GetShapeFct(double* Shape, double* m)
{ 
   GetDualMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
DualMortarQuad4::GetShapeFctVal(double* Shape, double* m)
{ 
   GetDualMortarShapeFct(Shape, m); 
}
