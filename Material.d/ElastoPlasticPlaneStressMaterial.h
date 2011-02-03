#ifndef ELASTOPLASTICPLANESTRESSMATERIAL
#define ELASTOPLASTICPLANESTRESSMATERIAL

#include <vector>

class ElastoPlasticPlaneStressMaterial
{
 public:
  
  //! Compute the elastoplastic constitutive response.
  //! Returns true if calculations went well and false otherwise.
  //! \param Fnp1 Input. Deformation gradient at new state of material. Size 9x1.
  //! \param CauchyStress Output. Has size 9x1.
  //! \param Cep Output. Algorithmic elastoplastic tangent. If requested, has size 81x1.
  //! \param UpdateFlag Input. Material state updated if true. Note that by default, material state is updated.
  virtual bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1, 
                                                        std::vector<double> * CauchyStress, 
                                                        std::vector<double> * Cep = 0, 
                                                        const bool UpdateFlag = true) = 0;
  
  //! Returns equivalent plastic strain in material
  virtual double GetMaterialEquivalentPlasticStrain() const = 0;

};

#endif
