#ifndef _SFEMINPC_H
#define _SFEMINPC_H

#include <Math.d/Vector.h>
#include <Utils.d/MyComplex.h>
#include <Sfem.d/Sfem.h>
#include <Driver.d/Domain.h>
#include<Math.d/matrix.h>

class Domain;
extern Domain *domain;
extern Sfem *sfem;
	
using namespace std;

template <class Scalar, class VecType>
class SfemInpc : public Sfem {
  int n;
//  double ** E_elem;
 public:
  SfemInpc() { P=sfem->getP(); L= sfem->getL(); ndim=sfem->getndim(); output_order = sfem->getoutput_order(); makealpha(); };
  ~SfemInpc() {};
  void setn(int n1) {n=n1;};
  void printCoefs(VecType* psi_u) {}; // print the coefficients in a file
  void computeMean(VecType* psi_u, VecType *mean_u);
  void computeStdDev(VecType* psi_u, VecType *sdev_u);
  void computePdf(int seed, VecType* psi_u, VecType *realz_u);
//  Scalar  computePdf(int seed, int dofno);
};


#ifdef _TEMPLATE_FIX_
#include <Sfem.d/SfemInpc.C>
#endif

#endif
