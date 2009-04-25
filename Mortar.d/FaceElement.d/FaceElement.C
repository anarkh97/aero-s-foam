// ------------------------------------------------------------
// HB -  08/11/03
// ------------------------------------------------------------
// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// FEM headers
#include <Utils.d/resize_array.h>

#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS
// -----------------------------------------------------------------------------------------------------
/*FaceElement::FaceElement():FFI(0) 
{ 
  Area = 0.0; 
  nFFI = 0;
}
*/

// -----------------------------------------------------------------------------------------------------
//                                           FFI METHODS 
// -----------------------------------------------------------------------------------------------------
/*
void
FaceElement::AddPtrFFI(FFIPolygon* _FFI)
{
  FFI[nFFI++] = _FFI;
}

int
FaceElement::nFFIs() { return nFFI; }

FFIPolygon*
FaceElement::GetPtrFFI(int i)
{
  return FFI[i];
}

void
FaceElement::printFFI()
{
  fprintf(stderr," -> nFFI = %d\n",nFFI);
}
*/

// -----------------------------------------------------------------------------------------------------
//                                        MAPPING & SHAPE FCT METHODS 
// -----------------------------------------------------------------------------------------------------
double* 
FaceElement::ViewRefCoords()
{
  fprintf(stderr," *** ERROR: method ViewRefCoords() NOT implemented for face el. type %d\n. Abort",
  GetFaceElemType());
  exit(1);
  return NULL;    //CRW
}

/*
double* 
FaceElement::ViewRefCoords()
{
 switch(GetFaceElemType())
 {
   case FaceElement::QUADFACEL4:
     return(&FaceQuad4::RefCoords);

   case FaceElement::QUADFACEQ8:
     return(&FaceQuad8::RefCoords);

   case FaceElement::QUADFACEQ9:
     return(&FaceQuad9::RefCoords);

   case FaceElement::QUADFACEC12:
     return(&FaceQuad12::RefCoords);

   case FaceElement::TRIFACEL3:
     return(&FaceTri3::RefCoords);

   case FaceElement::TRIFACEQ6:
     return(&FaceTri6::RefCoords);

   case FaceElement::TRIFACEC10:
     return(&FaceTri10::RefCoords);

   default:
     fprintf(stderr," *** ERROR: method ViewRefCoords() NOT implemented for face el. type %d\n. Abort",
     GetFaceElemType()));
     exit(-1);
     return;
  }
}
*/
