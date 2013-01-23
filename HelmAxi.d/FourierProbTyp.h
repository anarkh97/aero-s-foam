#ifndef _FOURIERPROBTYPE_H_
#define _FOURIERPROBTYPE_H_

class FourierStatic;
template<class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;

class FourierSolver {
     FourierStatic *probDesc;
     ComplexVector sol_n;
     ComplexVector rhs;
     ComplexVector *globalsol;

   public:

     FourierSolver();
     FourierSolver(FourierStatic *PrbD); 
     void solve();
};

#endif
