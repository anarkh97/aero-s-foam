#ifndef _NLMEMBRANE_H_
#define _NLMEMBRANE_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
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
  void getEandB(typename TwoDTensorTypes<n>::StrainTensor &e,
                typename TwoDTensorTypes<n>::BTensor &B,
                typename TwoDTensorTypes<n>::GradUTensor &gradU,
                typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
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
    void getEandB(typename TwoDTensorTypes<n>::StrainTensor &e,
                  typename TwoDTensorTypes<n>::BTensor &B,
                  typename TwoDTensorTypes<n>::GradUTensor &gradU,
                  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk);
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
    double interpolateScalar(double *_q, double _xi[3]);
};

class NLMembrane : public GenGaussIntgElement<TwoDTensorTypes<9> >
{
    int n[3];
    NLMaterial *material, *linearMaterial;
    double *cCoefs;
    double *cFrame;
    bool useDefaultMaterial;
    PressureBCond *pbc;

  protected:
    int getNumGaussPoints();
    void getGaussPointAndWeight(int i, double *point, double &weight);
    GenShapeFunction< TwoDTensorTypes<9> > *getShapeFunction();
    GenStrainEvaluator<TwoDTensorTypes<9> > *getGenStrainEvaluator();
    NLMaterial *getMaterial();
    NLMaterial *getLinearMaterial();
    void rotateConstitutiveMatrix2(CoordSet &cs, double C[6][6], double alpha[6]);

  public:
    NLMembrane(int *nd);
    ~NLMembrane();
    void setPreLoad(std::vector<double> &_preload) { preload = _preload; }
    std::vector<double> getPreLoad() { return preload; }
    int numNodes() { return 3; }
    int numDofs() { return 9; }
    void renum(int *);
    void renum(EleRenumMap&);
    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int* nodes(int * = 0);
    void updateStates(Node *nodes, double *states, double *un, double *unp) {}
    void setProp(StructProp *p, bool _myProp = false);
    void setCompositeData(int, int, double *, double *coefs, double *frame);
    double* setCompositeData2(int, int, double *, double *coefs, CoordSet &cs, double theta);
    void setMaterial(NLMaterial *);
    void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
    PressureBCond* getPressure() { return pbc; }
    void computePressureForce(CoordSet &cs, Vector& elPressureForce,
                              GeomState *gs=0, int cflg = 0, double t = 0);

    Corotator* getCorotator(CoordSet &, double *, int , int);
    int getTopNumber() { return 104; }
    FullSquareMatrix  stiffness(CoordSet& cs, double *k, int flg=1);
#ifdef USE_EIGEN3
    int getMassType() { return 2; } // both consistent and lumped
#else
    int getMassType() { return 0; } // lumped only
#endif
    FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
    double getMass(CoordSet&);

    void computeDisp(CoordSet &cs, State &state, const InterpPoint &, double *res, GeomState *gs);
    void getFlLoad(CoordSet &, const InterpPoint &,  double *flF, double *resF, GeomState *gs=0);

    PrioInfo examine(int sub, MultiFront *mf);
};

#include <Element.d/SuperElement.h>

class NLMembrane4 : public SuperElement
{
  public:
    NLMembrane4(int *nodenums);
    int  getTopNumber();
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    PrioInfo examine(int sub, MultiFront *mf);
};

#endif
