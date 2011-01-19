#ifndef ROM_BASIS_ORTHODRIVER_H
#define ROM_BASIS_ORTHODRIVER_H

class Domain;

class BasisOrthoDriver {
public:
  explicit BasisOrthoDriver(Domain *);

  void solve();

private:
  void preProcess();

  Domain *domain_;

  // Disallow copy and assignment
  BasisOrthoDriver(const BasisOrthoDriver &);
  BasisOrthoDriver &operator=(const BasisOrthoDriver &);
};

#endif /* ROM_BASIS_ORTHODRIVER_H */
