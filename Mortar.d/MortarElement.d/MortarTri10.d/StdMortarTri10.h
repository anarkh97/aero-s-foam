// ---------------------------------------------------------------- 
// HB - 05/24/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARTRI10_H_
#define _STDMORTARTRI10_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class FaceElement;

class StdMortarTri10: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarTri10();
        StdMortarTri10(FaceElement*);
        StdMortarTri10(double, FaceElement*);  
 
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
