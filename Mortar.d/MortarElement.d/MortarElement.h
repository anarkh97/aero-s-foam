// ---------------------------------------------------------------- 
// HB - 06/09/03
// HB - Modeified 02/28/04
// ---------------------------------------------------------------- 
#ifndef _MORTARELEMENT_H_ 
#define _MORTARELEMENT_H_ 

#include <cstdio>
//#include <iostream.h> 
//#include <Mortar.d/FaceElement.d/FaceElement.h>

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
};
#endif

