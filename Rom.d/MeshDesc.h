#ifndef ROM_MESHDESC_H
#define ROM_MESHDESC_H

#include "RenumberingUtils.h"

#include <Element.d/Element.h>
#include <Driver.d/StructProp.h>

#include <vector>
#include <map>
#include <iterator>
#include <cstring>
#include <utility>

class Domain;
class GeoSource;

class EFrameData;
class Attrib;
class NLMaterial;
class BCond;

namespace Rom {

class MeshDesc {
public:
  const CoordSet &nodes() const { return nodes_; }
  const Elemset &elements() const { return elements_; }

  const std::vector<EFrameData> &elemFrames() const { return elemFrames_; }
  const SPropContainer &properties() const { return *properties_; } 
  const std::vector<Attrib> &attributes() const { return attributes_; }
  const std::map<int, NLMaterial *> &materialLaws() const { return *materialLaws_; }
  const std::map<int, int> &materialLawMapping() const { return materialLawMapping_; }

  const std::vector<BCond> &dirichletBConds() const { return dirichletBConds_; }
  const std::vector<BCond> &neumannBConds() const { return neumannBConds_; }

  const std::vector<BCond> &initDisp() const { return initDisp_; }
  const std::vector<BCond> &initVel() const { return initVel_; }

  const std::map<int, double> &elemPressures() const { return elemPressures_; }
  
  // Only for a reduced mesh
  const std::vector<int> &sampleNodeIds() const { return sampleNodeIds_; }

  typedef std::map<int, Attrib> AttribContainer;

  // Creates a reduced mesh from the default one
  MeshDesc(Domain *, GeoSource *, const SampledMeshRenumbering &);

private:
  CoordSet nodes_; 
  Elemset elements_;

  std::vector<EFrameData> elemFrames_;
  std::vector<Attrib> attributes_;
  const SPropContainer *properties_;
  std::map<int, int> materialLawMapping_;
  const std::map<int, NLMaterial *> *materialLaws_;

  std::vector<BCond> dirichletBConds_;
  std::vector<BCond> neumannBConds_;
  std::vector<BCond> initDisp_; 
  std::vector<BCond> initVel_;
  std::map<int, double> elemPressures_;

  const std::vector<int> sampleNodeIds_;

  // Disallow copy & assignment
  MeshDesc(const MeshDesc &);
  MeshDesc &operator=(const MeshDesc &);
};

std::ostream &
operator<<(std::ostream &, const MeshDesc &);

} /* end namespace Rom */

#endif /* ROM_MESHDESC_H */
