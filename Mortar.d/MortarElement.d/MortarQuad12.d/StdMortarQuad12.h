// ---------------------------------------------------------------- 
// HB - 05/24/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD12_H_
#define _STDMORTARQUAD12_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

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

	// Shape fct methods
	// ~~~~~~~~~~~~~~~~~
	// -> local methods
	void GetStdMortarShapeFct(double* Shape, double* m);
	void GetShapeFct(double* Shape, double* m);

	// -> implementation of virtual methods
	void GetShapeFctVal(double* Shape, double* m);
};
#endif
