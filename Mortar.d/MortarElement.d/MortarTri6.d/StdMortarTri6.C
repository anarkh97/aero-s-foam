// ----------------------------------------------------------------
// HB - 08/25/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri6.d/StdMortarTri6.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarTri6::StdMortarTri6() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarTri6::StdMortarTri6(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarTri6::StdMortarTri6(double _area, FaceElement* FaceElem)
{
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}


// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarTri6::nNodes() { return(6); }

int
StdMortarTri6::nMortarShapeFct() { return(6); }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
StdMortarTri6::GetStdMortarShapeFct(double* Shape, double* m)
{
   GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

void
StdMortarTri6::GetShapeFct(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarTri6::GetShapeFctVal(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}
