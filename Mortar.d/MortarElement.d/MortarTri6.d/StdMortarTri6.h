// ---------------------------------------------------------------- 
// HB - 05/13/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARTRI6_H_
#define _STDMORTARTRI6_H_

//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class FaceElement;

class StdMortarTri6: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarTri6();
        StdMortarTri6(FaceElement*);
        StdMortarTri6(double, FaceElement*);  
 
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
