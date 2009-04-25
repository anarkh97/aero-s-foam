#ifndef _3DSHAPEFUNCTION_H_
#define _3DSHAPEFUNCTION_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Math.d/Vector.h>
class Node;

class ShapeFunction {
	
   protected:

     int numdofs;
     
   public:
     ShapeFunction(int n) : numdofs(n) {}
     virtual Tensor *getGradUInstance();
     virtual Tensor *getDgradUDqkInstance();
     virtual Tensor *getValInstance() = 0;
     virtual void getLocalDerivatives(Tensor *localDerivatives, double xi[3]) = 0;
     virtual void getValues(Tensor *gradU, Tensor *dgradUdqk, 
                  double *jac, double xi[3]) = 0;
     virtual void getGlobalGrads(Tensor *gradU, Tensor *dgradUdqk, double *jac, 
                  Node *nodes, double xi[3],Vector &disp);
     virtual void getGradU(Tensor *gradU, Node *nodes, double xi[3],Vector &disp);
};

class HexahedralShapeFunction : public ShapeFunction {

    public:
      HexahedralShapeFunction() : ShapeFunction(24) {}
      void getLocalDerivatives(Tensor *localDerivatives, double xi[3]);
      void getValues(Tensor *gradU, Tensor *dgradUdqk, double *jac, double xi[3]) {}
      Tensor *getValInstance(){return 0;}
};


template <class TensorTypes>
class GenShapeFunction {
   public:
     GenShapeFunction() {}
     virtual void getGlobalGrads(typename TensorTypes::GradUTensor &gradU,
      typename TensorTypes::GradUDerivTensor &dgradUdqk, double *jac,
                  Node *nodes, double xi[3],Vector &disp) = 0;
     virtual void getGradU(typename TensorTypes::GradUTensor &gradU,
          Node *nodes, double xi[3],Vector &disp) = 0;
};

#endif
