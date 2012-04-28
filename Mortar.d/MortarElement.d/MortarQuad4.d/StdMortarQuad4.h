// ---------------------------------------------------------------- 
// HB - 06/09/03
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD4_H_
#define _STDMORTARQUAD4_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class FaceElement;

class StdMortarQuad4: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarQuad4();
        StdMortarQuad4(FaceElement*);
        StdMortarQuad4(double, FaceElement*);  

        virtual ~StdMortarQuad4() { }
 
	// Get methods
	// ~~~~~~~~~~~
	// -> implementation of virtual methods
	int nNodes();
	int nMortarShapeFct();

	// Shape fct methods
	// ~~~~~~~~~~~~~~~~~
	// -> local methods
        template<typename Scalar>
	  void GetStdMortarShapeFct(Scalar* Shape, Scalar* m);
	void GetShapeFct(double* Shape, double* m);

	// -> implementation of virtual methods
	void GetShapeFctVal(double* Shape, double* m);
#if (MAX_MORTAR_DERIVATIVES > 0)
        void GetShapeFctVal(MadDouble* Shape, MadDouble* m);
#endif
};

template<typename Scalar>
void
StdMortarQuad4::GetStdMortarShapeFct(Scalar* Shape, Scalar* m)
{
   GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

#endif
