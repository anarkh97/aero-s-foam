#ifndef PITA_HALFSLICENETWORKTOPOLOGY_H
#define PITA_HALFSLICENETWORKTOPOLOGY_H

#include <vector>
#include <map>

#include "Fwk.h"
#include "Types.h"

namespace Pita {

class HalfSliceNetworkTopology : public Fwk::PtrInterface<HalfSliceNetworkTopology> {
public:
  typedef Fwk::Ptr<HalfSliceNetworkTopology> Ptr;
  typedef Fwk::Ptr<const HalfSliceNetworkTopology> PtrConst;
  
  template <typename T>
  class GenIteratorConst {
  public:
    const T & operator*() const { return *current_; }
    const T * operator->() const { return current_.operator->(); }
    GenIteratorConst<T> & operator++() { ++current_; return *this; }
    GenIteratorConst<T> operator++(int) { GenIteratorConst<T> tmp(*this); (*this)++; return tmp; }
    operator bool() const { return current_ != end_; }

    // Default SeedIteratorConst(const SliceIteratorConst &);
    // Default SeedIteratorConst & operator=(const SliceIteratorConst &);
  private:
    friend class HalfSliceNetworkTopology;
    
    typedef typename std::vector<T>::const_iterator ItConst;
    ItConst current_, end_;

    explicit GenIteratorConst(ItConst begin, ItConst end) :
      current_(begin), end_(end)
    {}
  };

  typedef GenIteratorConst<Hs::SliceId> SliceIteratorConst;
  typedef GenIteratorConst<Hs::SeedId> SeedIteratorConst;

  SliceIteratorConst writer(const Hs::SeedId & id) const;
  SliceIteratorConst reader(const Hs::SeedId & id) const;
  SeedIteratorConst seed(const Hs::SliceId & id) const;
  
  static HalfSliceNetworkTopology::Ptr New() {
    return new HalfSliceNetworkTopology();
  }

protected:
  friend class GenIteratorConst<Hs::SeedId>;
  friend class GenIteratorConst<Hs::SliceId>;

  HalfSliceNetworkTopology() {}

private:
  typedef std::vector<Hs::SliceId> SliceContainer;
  typedef std::vector<Hs::SeedId> SeedContainer;

  mutable std::map<Hs::SeedId, SliceContainer> writer_, reader_;
  mutable std::map<Hs::SliceId, SeedContainer> seed_;

  template <typename K, typename T>
    GenIteratorConst<T> createIterator(std::map<K, std::vector<T> > & m, const K & id,
                                       void (* const buildFunc)(const K &, std::vector<T> &)) const;
};

} // end namespace Pita

#endif /* PITA_HALFSLICENETWORKTOPOLOGY_H */
