#ifndef ROM_DRIVERINTERFACE_H
#define ROM_DRIVERINTERFACE_H

class Domain;

class RomDriverInterface {
public:
  virtual void solve() = 0;
  virtual ~RomDriverInterface() {}

protected:
  RomDriverInterface() {}

private:
  // Disallow copy & assignment
  RomDriverInterface(const RomDriverInterface &);
  RomDriverInterface &operator=(const RomDriverInterface &);
};

// Concrete class instantiation
extern RomDriverInterface *basisOrthoDriverNew(Domain *);
extern RomDriverInterface *meshSamplingDriverNew(Domain *);

#endif /* ROM_DRIVERINTERFACE_H */
