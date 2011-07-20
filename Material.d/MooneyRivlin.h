#include <vector>
#include <string>
#include <cmath>
#include "Material.h"


class MooneyRivlin: public SimpleMaterial
{
 public:
  MooneyRivlin(double Mu1Input,double Mu2Input, double KappaInput, double rhoInput)
    :SimpleMaterial(rhoInput), Mu1(Mu1Input), Mu2(Mu2Input), Kappa(KappaInput) {}
  virtual ~MooneyRivlin() {}
  MooneyRivlin(const MooneyRivlin &NewMat): 
    SimpleMaterial(NewMat), Mu1(NewMat.Mu1), Mu2(NewMat.Mu2), Kappa(NewMat.Kappa) {}
  virtual MooneyRivlin * Clone() const { return new MooneyRivlin(*this); }
  
  bool GetConstitutiveResponse(const std::vector<double> * strain,
			       std::vector<double> * stress,
			       std::vector<double> * tangents = 0) const;

  const std::string GetMaterialName() const { return "MooneyRivlin"; }

  double GetSoundSpeed(void) const { return sqrt(2.*Kappa-4./3.*(2.0*(Mu1+Mu2)) + 2.0 * (2.0*(Mu1+Mu2))/GetDensityInReference());}; 

 private:
  // Lame constants
  double Mu1;
  double Mu2;
  double Kappa;
};