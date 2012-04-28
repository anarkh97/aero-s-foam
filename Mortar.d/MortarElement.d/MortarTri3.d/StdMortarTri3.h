// ---------------------------------------------------------------- 
// HB - 08/25/03
// ---------------------------------------------------------------- 
#ifndef _STDMORTARTRI3_H_
#define _STDMORTARTRI3_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class FaceElement;

class StdMortarTri3: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarTri3();
        StdMortarTri3(FaceElement*);
        StdMortarTri3(double, FaceElement*);  

        virtual ~StdMortarTri3() { }
 
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
StdMortarTri3::GetStdMortarShapeFct(Scalar* Shape, Scalar* m)
{
   GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

#endif
