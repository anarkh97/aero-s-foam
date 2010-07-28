#ifndef _AUXDEFS_H_
#define _AUXDEFS_H_

#include<Element.d/Element.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/SimpleMat.h>
#include <Element.d/NonLinearity.d/HypoElasticMat.h>
#include <Element.d/NonLinearity.d/ElastoViscoPlasticMat.h>
#include <Element.d/NonLinearity.d/J2PlasticityMat.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>

#ifndef TFLOP
#include <HelmAxi.d/FourierHelmBCs.h>
#endif

#include <HelmAxi.d/ScatterData.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <HelmAxi.d/MPCData.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <HelmAxi.d/Line2AxiSommer.h>

#include <Parser.d/DecInit.h>

class Domain;
class GeoSource;
class Sfem;
template <class Scalar, class VecType> class SfemNonInpc;
template <class Scalar, class VecType> class SfemInpc;

typedef struct {
   int num;
   double xyz[3];
} NumedNode;
   
typedef struct {
   int num;
   int nd[125];
} NumList;

struct NumedElem {
  int num;
  Element *elem;
};

class BCList {
  public:
    BCond *d;
    int n;
    int maxbc;

    BCList();
    void add(BCond &);
    void add(int nd, int dof, double val) 
      { BCond bc; bc.nnum = nd; bc.dofnum =dof; bc.val = val; add(bc); }
};

class ComplexBCList {
  public:
    ComplexBCond *d;
    int n;
    int maxbc;

    ComplexBCList();
    void add(ComplexBCond &);
};

struct FrameData {
  int num;
  double d[9];
};

struct LayerData {
  int lnum;
  int matid;
  double d[9];
};

struct ConstraintOptions {
  bool lagrangeMult;
  double penalty;
};

extern double fetiHTol;
extern int fetiHIterMax;

extern DecInit * decInit;

extern Domain *domain;
extern GeoSource *geoSource;

extern int line_num;
extern int numColumns;

int yylex(void);

void yyerror(const char*);

// dynamic weight for every type of elements
// element numbers are those assigned in AddEle.C
extern map<int, double > weightList;  

#endif
