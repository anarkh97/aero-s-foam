// ----------------------------------------------------------------
// HB - 05/24/05
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri10.d/StdMortarTri10.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarTri10::StdMortarTri10() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarTri10::StdMortarTri10(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarTri10::StdMortarTri10(double _area, FaceElement* FaceElem)
{
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}


// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarTri10::nNodes() { return(10); }

int
StdMortarTri10::nMortarShapeFct() { return(10); }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
StdMortarTri10::GetStdMortarShapeFct(double* Shape, double* m)
{
   GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

void
StdMortarTri10::GetShapeFct(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarTri10::GetShapeFctVal(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}
