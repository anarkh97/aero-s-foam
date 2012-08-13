// ----------------------------------------------------------------
// HB - 06/09/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarQuad4.d/StdMortarQuad4.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarQuad4::StdMortarQuad4() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarQuad4::StdMortarQuad4(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarQuad4::StdMortarQuad4(double area_, FaceElement* FaceElem)
{
  SetArea(area_);
  SetPtrMasterFace(FaceElem);
}


// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarQuad4::nNodes() { return 4; }

int
StdMortarQuad4::nMortarShapeFct() { return 4; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
StdMortarQuad4::GetShapeFct(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarQuad4::GetShapeFctVal(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

#if (MAX_MORTAR_DERIVATIVES > 0)
void
StdMortarQuad4::GetShapeFctVal(ActiveDouble* Shape, ActiveDouble* m)
{
   GetStdMortarShapeFct(Shape, m);
}
#endif
