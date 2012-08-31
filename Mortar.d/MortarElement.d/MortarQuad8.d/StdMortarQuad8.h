// ---------------------------------------------------------------- 
// HB - 05/13/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD8_H_
#define _STDMORTARQUAD8_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class FaceElement;

class StdMortarQuad8: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarQuad8();
        StdMortarQuad8(FaceElement*);
        StdMortarQuad8(double, FaceElement*);  
 
	// Get methods
	// ~~~~~~~~~~~
	// -> implementation of virtual methods
	int nNodes();
	int nMortarShapeFct();

	// Shape fct methods
	// ~~~~~~~~~~~~~~~~~
	// -> local methods
	void GetStdMortarShapeFct(double* Shape, double* m);
	void GetShapeFct(double* Shape, double* m);

	// -> implementation of virtual methods
	void GetShapeFctVal(double* Shape, double* m);
};
#endif
