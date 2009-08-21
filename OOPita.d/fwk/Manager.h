#ifndef FWK_MANAGER_H
#define FWK_MANAGER_H

#include "NamedInterface.h"
#include "Exception.h"
#include "Notifiee.h"
#include <map>

namespace Fwk {

/* Generic Manager Implementation
** T => Managed type
** K => Key                       */
template <typename T, typename K>
class GenManagerImpl {
public:
  typedef std::map<K, Ptr<T> > InstanceMap;
  typedef typename InstanceMap::size_type InstanceCount;
  typedef typename InstanceMap::const_iterator IteratorConst;
  typedef typename InstanceMap::iterator Iterator;
  typedef typename InstanceMap::value_type Pair;

  /*class IteratorConst {
  public:
    const Pair & operator*() const { return *it_; }
    const Pair * operator->() const { return it_.operator->(); }
    IteratorConst & operator++() { ++it_; return *this; }
    IteratorConst operator++(int) { IteratorConst tmp(*this); ++(*this); return tmp; }

  private:
    typename InstanceMap::const_iterator it_;
  };*/

  T * instance(const K & key) const;
  InstanceCount instanceCount() const { return instance_.size(); }
  
  T * instanceNew(const K & key);
  void instanceDel(const K & key) { instance_.erase(key); } 
  
  IteratorConst instanceBegin() const { return instance_.begin(); }
  IteratorConst instanceEnd() const { return instance_.end(); }

  Iterator instanceBegin() { return instance_.begin(); }
  Iterator instanceEnd() { return instance_.end(); }

protected:
  virtual ~GenManagerImpl() {}
  virtual T * createNewInstance(const K & key) = 0;

private:
  InstanceMap instance_;
};

template <typename T, typename K>
T *
GenManagerImpl<T, K>::instance(const K & key) const {
  typename InstanceMap::const_iterator it = instance_.find(key);
  return (it != instance_.end()) ? it->second.ptr() : NULL;
}

template <typename T, typename K>
T *
GenManagerImpl<T, K>::instanceNew(const K & key) {
  // Find insertion point
  typename InstanceMap::iterator it = instance_.lower_bound(key);
  if (it != instance_.end() && it->first == key)
    throw NameInUseException();
 
  // Build and add new instance
  Ptr<T> newInstance = createNewInstance(key);
  instance_.insert(it, std::make_pair(key, newInstance));

  return newInstance.ptr();
}

typedef GenManagerImpl<NamedInterface, String> NamedInterfaceManagerImpl;

template <typename T>
class GenNamedInterfaceManager : public PtrInterface<GenNamedInterfaceManager<T> >, private NamedInterfaceManagerImpl {
public:
  typedef Fwk::Ptr<GenNamedInterfaceManager<T> > Ptr;
  typedef Fwk::Ptr<const GenNamedInterfaceManager<T> > PtrConst;
 
  typedef NamedInterfaceManagerImpl::InstanceCount InstanceCount;

  T * instance(const String & name) const { return static_cast<T *>(NamedInterfaceManagerImpl::instance(name)); } 
  InstanceCount instanceCount() const { return NamedInterfaceManagerImpl::instanceCount(); }
  
  T * instanceNew(const String & name) { return static_cast<T *>(NamedInterfaceManagerImpl::instanceNew(name)); }
  void instanceDel(const String & name) { NamedInterfaceManagerImpl::instanceDel(name); }

protected:
  GenNamedInterfaceManager() {}
  GenNamedInterfaceManager(const GenNamedInterfaceManager<T> &); // No implementation
  GenNamedInterfaceManager<T> & operator=(const GenNamedInterfaceManager<T> &); // No implementation

  virtual T * createNewInstance(const String & key) = 0;
};

typedef GenNamedInterfaceManager<NamedInterface> NamedInterfaceManager;

/* Generic Manager Interface */
template <typename PtrT, typename K>
class GenManagerInterface : public Fwk::PtrInterface<GenManagerInterface<PtrT, K> > {
public:
  typedef Fwk::Ptr<GenManagerInterface<PtrT, K> > Ptr;
  typedef Fwk::Ptr<const GenManagerInterface<PtrT, K> > PtrConst;

  // Read accessors
  virtual PtrT instance(const K & key) const = 0;
  virtual size_t instanceCount() const = 0;

  // Mutators
  virtual PtrT instanceNew(const K & key) = 0;
  virtual void instanceDel(const K & key) = 0;

  // Notification support
  /*class NotifieeConst : public BaseNotifiee<const GenManagerInterface, NotifieeConst> {
  public:
    typedef Fwk::Ptr<NotifieeConst> Ptr;
    typedef Fwk::Ptr<const NotifieeConst> PtrConst;

    virtual void onInstanceNew(const K & key) {}
    virtual void onInstanceDel(const K & key) {}

  protected:
    explicit NotifieeConst(const GenManagerInterface<PtrT, K> * notifier) :
      BaseNotifiee<const GenManagerInterface<PtrT, K> >(notifier)
    {}
  };

  NotifieeConst * lastNotifiee() const { return notifiee_; }
  virtual void lastNotifieeIs(NotifieeConst * n) { setNotifiee(n); }

protected:
  GenManagerInterface() : notifiee_() {}

  void setNotifiee(NotifieeConst * n) { notifiee_ = n; }

private:
  NotifieeConst * notifiee_;*/
};

} // end namespace Fwk

#endif /* Fwk_MANAGER_H */
