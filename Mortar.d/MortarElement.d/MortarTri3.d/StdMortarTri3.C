// ----------------------------------------------------------------
// HB - 08/25/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri3.d/StdMortarTri3.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarTri3::StdMortarTri3() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarTri3::StdMortarTri3(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarTri3::StdMortarTri3(double _area, FaceElement* FaceElem)
{
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}


// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarTri3::nNodes() { return 3; }

int
StdMortarTri3::nMortarShapeFct() { return 3; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
StdMortarTri3::GetShapeFct(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarTri3::GetShapeFctVal(double* Shape, double* m)
{ 
   GetStdMortarShapeFct(Shape, m); 
}

#if (MAX_MORTAR_DERIVATIVES > 0)
void
StdMortarTri3::GetShapeFctVal(MadDouble* Shape, MadDouble* m)
{
   GetStdMortarShapeFct(Shape, m);
}
#endif

