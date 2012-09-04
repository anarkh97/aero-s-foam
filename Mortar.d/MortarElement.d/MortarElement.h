// ---------------------------------------------------------------- 
// HB - 06/09/03
// HB - Modeified 02/28/04
// ---------------------------------------------------------------- 
#ifndef _MORTARELEMENT_H_ 
#define _MORTARELEMENT_H_ 

#include <iostream>
//#include <iostream.h> 
//#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarAutoDiff.h>

class FaceElement;

class MortarElement{
  private:
	double Area;

        FaceElement* MasterFace; // ptr to supporting face element 

  public:
	// Constructors
	// ~~~~~~~~~~~~
	MortarElement(); 
	MortarElement(FaceElement*);
        MortarElement(double, FaceElement*); 
	
	// Destructors
	// ~~~~~~~~~~~
	virtual ~MortarElement(); 
 
        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

	// Set methods
	// ~~~~~~~~~~~
	// -> local
	void SetArea(double);
	void SetPtrMasterFace(FaceElement*);

	// Get methods
	// ~~~~~~~~~~~
	// -> local
	int nMasterFaceNodes();   
	FaceElement* GetPtrMasterFace(); 
	
	// -> interface 
	virtual int nNodes() ;
	virtual int nMortarShapeFct();

	// Shape fct methods
        // ~~~~~~~~~~~~~~~~~
	// -> interface 
	virtual void GetShapeFctVal(double* Shape, double* m);
#if (MAX_MORTAR_DERIVATIVES > 0)
        virtual void GetShapeFctVal(ActiveDouble* Shape, ActiveDouble* m)
          { std::cerr << "MortarElement::GetShapeFctVal(ActiveDouble* Shape, ActiveDouble* m) is not implemented\n"; }
#endif
};
#endif

