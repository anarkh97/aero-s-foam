#include <Corotational.d/Corotator.h>
#include <iostream>

void 
Corotator::formGeometricStiffness(GeomState&, CoordSet&, FullSquareMatrix&, double*) 
{ 
  std::cerr << " *** WARNING: Corotator::formGeometricStiffness(GeomState&, CoordSet&, FullSquareMatrix&, double*) is not implemented\n";
}

double* 
Corotator::getOriginalStiffness() 
{ 
  std::cerr << " *** WARNING: Corotator::getOriginalStiffness() is not implemented\n";
  return 0; 
}

void
Corotator::extractDeformations(GeomState&, CoordSet&, double*, int&) 
{
  std::cerr << " *** WARNING: Corotator::extractDeformations(GeomState&, CoordSet&, double*, int&) is not implemented\n"; 
}

void
Corotator::extractDeformations(GeomState&, CoordSet&, DComplex*, int&)
{
  std::cerr << " *** WARNING: Corotator::extractDeformations(GeomState&, CoordSet&, DComplex*, int&) is not implemented\n";
}

void
Corotator::getNLVonMises(Vector&, Vector&, GeomState&, CoordSet&, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMises(Vector&, Vector&, GeomState&, CoordSet&, int) is not implemented\n";
}

void
Corotator::getNLVonMises(ComplexVector&, Vector&, GeomState&, CoordSet&, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMises(ComplexVector&, Vector&, GeomState&, CoordSet&, int) is not implemented\n";
}

void
Corotator::getNLAllStress(FullM&, Vector&, GeomState&, CoordSet&, int)
{
  std::cerr << " *** WARNING: Corotator::getNLAllStress(FullM&, Vector&, GeomState&, CoordSet&, int) is not implemented\n";
}

double
Corotator::getElementEnergy(GeomState&, CoordSet&)
{
  std::cerr << " *** WARNING: Corotator::getElementEnergy(GeomState&, CoordSet&) is not implemented\n";
  return 0.0;
}

void
Corotator::extractRigidBodyMotion(GeomState&, CoordSet&, double*)
{
  std::cerr << " *** WARNING: Corotator::extractRigidBodyMotion(GeomState&, CoordSet&, double*) is not implemented\n";
}

