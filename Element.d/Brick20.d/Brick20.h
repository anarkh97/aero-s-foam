#ifndef _TWENTY_BRICK_H_
#define _TWENTY_BRICK_H_

#include <Element.d/Element.h>

#include <complex>
using std::complex;

class Brick20: public Element
{
    int nn[20];
    double *cCoefs;
    double *cFrame;
    NLMaterial *mat;

  public:
    Brick20(int*);
    ~Brick20();

    Element *clone();

    void renum(int *);
    void renum(EleRenumMap&);

    FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
    void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega);
    double getMass(CoordSet& cs);

    void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs);
    void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs=0);

    void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                     int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                      int surface=0, double *ndTemps=0);

    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int numDofs();

    int numNodes();
    int* nodes(int * = 0);

    int getTopNumber();

    PrioInfo examine(int sub, MultiFront *);
    int nDecFaces() { return 6; }
    int getDecFace(int iFace, int *fn);

    int getFace(int iFace, int *fn);

    void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame)
      { cCoefs = coefs; cFrame = frame; }

    double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
      { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                "              for Hexa20 el.\n"); return (double *) 0;
      }

    void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                          int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

    void setMaterial(NLMaterial *);
    int numStates();
    void initStates(double *st);
    Corotator *getCorotator(CoordSet &cs, double *kel, int=2, int=2);
};

#endif
