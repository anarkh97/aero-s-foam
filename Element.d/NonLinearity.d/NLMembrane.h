#ifndef _NLMEMBRANE_H_
#define _NLMEMBRANE_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainDispEvaluator.h>
#include <Math.d/TTensor.h>

template <int n>
struct TwoDTensorTypes
{
   static const int ndofs=n;
   typedef Stress2D StressTensor;
   typedef Stress2D StrainTensor;
   typedef SimpleTensor<StrainTensor, n> BTensor;
   typedef SimpleTensor<SimpleTensor<double,3>,2> GradUTensor;
   typedef SimpleTensor<GradUTensor, n> GradUDerivTensor;
   typedef SimpleTensor<BTensor,n> DBTensor;
   // Material tangent matrix tensor type
   typedef SymTensor<StressTensor,2> DTensor;
};

// Double contraction operator
template <int n>
SimpleTensor<Stress2D, n>
operator||(SimpleTensor<SimpleTensor<Stress2D,n>,n> &, SimpleTensor<Stress2D, n> &);

template <int n>
class LinearStrain2D : public GenStrainEvaluator<TwoDTensorTypes<n> >
{
  void getE(typename TwoDTensorTypes<n>::StrainTensor &e,
            typename TwoDTensorTypes<n>::GradUTensor &gradU);
  void getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                  typename TwoDTensorTypes<n>::BTensor &B, 
                  typename TwoDTensorTypes<n>::DBTensor &DB,
                  typename TwoDTensorTypes<n>::GradUTensor &gradU, 
                  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
};

template <int n>
class GLStrain2D : public GenStrainEvaluator<TwoDTensorTypes<n> > 
{
  public:
    void getE(typename TwoDTensorTypes<n>::StrainTensor &e,
              typename TwoDTensorTypes<n>::GradUTensor &gradU);
    void getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                    typename TwoDTensorTypes<n>::BTensor &B, 
                    typename TwoDTensorTypes<n>::DBTensor &DB,
                    typename TwoDTensorTypes<n>::GradUTensor &gradU, 
                    typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
    bool isNonLinear() { return true; }
};

class TriMembraneShapeFunct : public GenShapeFunction< TwoDTensorTypes<9> >
{
  public:
    TriMembraneShapeFunct() : GenShapeFunction< TwoDTensorTypes<9> >() {}
    void getGlobalGrads(Grad2D &gradU, Grad2DDeriv9 &dGradUdqk, double *_jac,
                        Node *nodes, double xi[3], Vector &disp);
    void getGradU(Grad2D &gradU,
                  Node *nodes, double xi[3], Vector &disp);
};

class NLMembrane : public GenGaussIntgElement<TwoDTensorTypes<9> >
{
    int n[3];
    NLMaterial *material;
    bool linearKinematics;

  protected:
    int getNumGaussPoints();
    void getGaussPointAndWeight(int i, double *point, double &weight);
    GenShapeFunction< TwoDTensorTypes<9> > *getShapeFunction();
    StrainEvaluator *getStrainEvaluator();
    GenStrainEvaluator<TwoDTensorTypes<9> > *getGenStrainEvaluator();
    NLMaterial *getMaterial();

  public:
    NLMembrane(int *nd, bool isLinKin);
    int numNodes() { return 3; }
    int numDofs() { return 9; }
    void renum(int *);
    void   markDofs(DofSetArray &);
    int*   dofs(DofSetArray &, int *p=0);
    int*   nodes(int * = 0);
    void updateStates(Node *nodes, double *states, double *un, double *unp) {}
    void setProp(StructProp *);
    void setMaterial(NLMaterial *);
    void computePressureForce(Node *cs,Vector& elPressureForce,
                              double *gs=0, int cflg = 0);
};

#endif
