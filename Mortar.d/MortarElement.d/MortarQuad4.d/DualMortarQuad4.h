// ---------------------------------------------------------------- 
// HB - 08/13/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARQUAD4_H_
#define _DUALMORTARQUAD4_H_

//#include <Math.d/matrix.h>
//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarQuad4: public MortarElement {
  private:
	FullM* Alpha; // coeffs for computing dual shape fcts 
                      // from the std ones: Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarQuad4();
        DualMortarQuad4(FaceElement*);
        DualMortarQuad4(FaceElement*, CoordSet&);
        
	DualMortarQuad4(double, FaceElement*);  
        DualMortarQuad4(double, FaceElement*, CoordSet&);  
        
	// Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarQuad4();

        // Initialization & clean/clear methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();
 
	// Get methods
	// ~~~~~~~~~~~
	// -> implementation of virtual methods
	int nNodes();
	int nMortarShapeFct();

	// Shape fct methods
	// ~~~~~~~~~~~~~~~~~
	// -> local methods
        void ComputeDualCoeffs(CoordSet&);
	void GetDualMortarShapeFct(double* Shape, double* m);
	void GetShapeFct(double* Shape, double* m);

	// -> implementation of virtual methods
	void GetShapeFctVal(double* Shape, double* m);
};
#endif
