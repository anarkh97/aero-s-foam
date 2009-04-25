// ----------------------------------------------------------------
// HB - 06/09/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarQuad8.d/StdMortarQuad8.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarQuad8::StdMortarQuad8() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarQuad8::StdMortarQuad8(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarQuad8::StdMortarQuad8(double area_, FaceElement* FaceElem)
{
  SetArea(area_);
  SetPtrMasterFace(FaceElem);
}


// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarQuad8::nNodes() { return(8); }

int
StdMortarQuad8::nMortarShapeFct() { return(8); }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
StdMortarQuad8::GetStdMortarShapeFct(double* Shape, double* m)
{
   GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

void
StdMortarQuad8::GetShapeFct(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarQuad8::GetShapeFctVal(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}
