// ---------------------------------------------------------------- 
// HB - 05/24/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD12_H_
#define _STDMORTARQUAD12_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.h>

class FaceElement;

class StdMortarQuad12: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarQuad12();
        StdMortarQuad12(FaceElement*);
        StdMortarQuad12(double, FaceElement*);  

        // Get methods
        // ~~~~~~~~~~~
        // -> implementation of virtual methods
        int nNodes();
        int nMortarShapeFct();
        bool GetDualFlag() { return false; }

        // Shape fct methods
        // ~~~~~~~~~~~~~~~~~
        // -> local methods
        template<typename Scalar>
          void GetShapeFctVal(Scalar* Shape, Scalar* m);
        template<typename Scalar>
          void GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar>
          void Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m);

        // -> implementation of virtual methods
        void GetShapeFctVal(double* Shape, double* m);
        void GetdShapeFct(double* dShapex, double* dShapey, double* m);
        void Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m);
};

template<typename Scalar>
void
StdMortarQuad12::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  static_cast<FaceQuad12*>(GetPtrMasterFace())->GetShapeFctVal(Shape, m);
}

template<typename Scalar>
void
StdMortarQuad12::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  static_cast<FaceQuad12*>(GetPtrMasterFace())->GetdShapeFct(dShapex, dShapey, m);
}

template<typename Scalar>
void
StdMortarQuad12::Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m)
{
  static_cast<FaceQuad12*>(GetPtrMasterFace())->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}

#endif
