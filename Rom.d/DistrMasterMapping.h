#ifndef ROM_DISTRMASTERMAPPING_H
#define ROM_DISTRMASTERMAPPING_H

#include "MasterMapping.h"
#include "DistrDomainUtils.h"

#include <Driver.d/SubDomain.h>

#include <vector>

namespace Rom {

class DistrMasterMapping {
public:
  int subDomainCount() const { return subMapping_.size(); }
  int localNodeCount() const { return localNodes_.size(); }
  int masterNodeCount() const { return masterNodes_.size(); }

  typedef std::vector<MasterMapping>::const_iterator SubMasterMappingIt;
  SubMasterMappingIt begin() const { return subMapping_.begin(); }
  SubMasterMappingIt end()   const { return subMapping_.end();   }

  typedef std::vector<int>::const_iterator NodeIndexIt;
  NodeIndexIt localNodeBegin() const { return localNodes_.begin(); }
  NodeIndexIt localNodeEnd()   const { return localNodes_.end();   }
  NodeIndexIt masterNodeBegin() const { return masterNodes_.begin(); }
  NodeIndexIt masterNodeEnd()   const { return masterNodes_.end();   }

  DistrMasterMapping() {}
  template <typename SubDomFwdIt>
  DistrMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast);

protected:
  std::vector<MasterMapping> subMapping_;
  std::vector<int> localNodes_;
  std::vector<int> masterNodes_;
};

template <typename SubDomFwdIt>
DistrMasterMapping::DistrMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast) {
  for (SubDomFwdIt subDomIt = subDomFirst; subDomIt != subDomLast; ++subDomIt) {
    SubDomain &sd = *subDomIt;
    const int * const globalBegin = sd.getGlNodes();
    const int * const globalEnd = globalBegin + sd.numNodes();

    localNodes_.insert(localNodes_.end(), globalBegin, globalEnd);

    std::vector<bool> masterNodes;
    master_node_flags(sd, std::back_inserter(masterNodes));
    subMapping_.push_back(MasterMapping(globalBegin, globalEnd, masterNodes.begin()));

    std::vector<bool>::const_iterator flagIt = masterNodes.begin();
    for (const int *globalIt = globalBegin; globalIt != globalEnd; ++globalIt) {
      if (*flagIt++) {
        masterNodes_.push_back(*globalIt);
      }
    }
  }
}

class DistrMpcMasterMapping : public DistrMasterMapping {
public:
  template <typename SubDomFwdIt>
  DistrMpcMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast);
};

template <typename SubDomFwdIt>
DistrMpcMasterMapping::DistrMpcMasterMapping(SubDomFwdIt subDomFirst, SubDomFwdIt subDomLast) {
  for (SubDomFwdIt subDomIt = subDomFirst; subDomIt != subDomLast; ++subDomIt) {
    SubDomain &sd = *subDomIt;
    const int * const globalBegin = sd.getGlMPCs();
    const int * const globalEnd = globalBegin + sd.getNumMpc();

    localNodes_.insert(localNodes_.end(), globalBegin, globalEnd);

    std::vector<bool> masterNodes(sd.getMpcMaster(), sd.getMpcMaster()+sd.getNumMpc());
    subMapping_.push_back(MasterMapping(globalBegin, globalEnd, masterNodes.begin()));

    std::vector<bool>::const_iterator flagIt = masterNodes.begin();
    for (const int *globalIt = globalBegin; globalIt != globalEnd; ++globalIt) {
      if (*flagIt++) {
        masterNodes_.push_back(*globalIt);
      }
    }
  }
}

} // end namespace Rom


#endif /* ROM_DISTRMASTERMAPPING_H */
