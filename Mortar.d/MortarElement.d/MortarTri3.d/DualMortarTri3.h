// ---------------------------------------------------------------- 
// HB - 08/26/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARTRI3_H_
#define _DUALMORTARTRI3_H_

//#include <Math.d/matrix.h>
//#include <Element.d/Element.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarTri3: public MortarElement {
  private:
        // !!! NO NEED TO COMPUTE THE Alpha COEFFS:
        // -> THEY ARE A PRIORI INDEPENDANT OF THE FACE ELEMENT GEOMETRY
	//double Alpha[3]; // coeffs for computing dual shape fcts 
                           // from the std ones: 
                           //  Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarTri3();
        DualMortarTri3(FaceElement*);
        DualMortarTri3(FaceElement*, CoordSet&);
        
	DualMortarTri3(double, FaceElement*);  
        DualMortarTri3(double, FaceElement*, CoordSet&);  
        
	// Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarTri3();
	
	// Set methods
	// ~~~~~~~~~~~
	// -> local methods
	//void SetDualCoeffs();
 
	// Get methods
	// ~~~~~~~~~~~
	// -> implementation of virtual methods
	int nNodes();
	int nMortarShapeFct();

	// Shape fct methods
	// ~~~~~~~~~~~~~~~~~
	// -> local methods
        //void ComputeDualCoeffs(CoordSet&);
	void GetDualMortarShapeFct(double* Shape, double* m);
	void GetShapeFct(double* Shape, double* m);

	// -> implementation of virtual methods
	void GetShapeFctVal(double* Shape, double* m);
};
#endif
