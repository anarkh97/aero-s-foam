#ifndef ROM_DRIVERINTERFACE_H
#define ROM_DRIVERINTERFACE_H

class Domain;

namespace Rom {

class DriverInterface {
public:
  virtual void solve() = 0;
  virtual ~DriverInterface() {}

protected:
  DriverInterface() {}

private:
  // Disallow copy & assignment
  DriverInterface(const DriverInterface &);
  DriverInterface &operator=(const DriverInterface &);
};

} /* end namespace Rom */

// Concrete class instantiation
extern Rom::DriverInterface *basisOrthoDriverNew(Domain *);
extern Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *);
extern Rom::DriverInterface *meshSamplingDriverNew(Domain *);

#endif /* ROM_DRIVERINTERFACE_H */
