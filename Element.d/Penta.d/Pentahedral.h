#ifndef _PENTAHEDRAL_H_
#define _PENTAHEDRAL_H_

#include <Element.d/Element.h>

class Pentahedral: public Element
{
    int nn[6];
    double *cCoefs;
    double *cFrame;
    NLMaterial *mat;

  public:
    Pentahedral(int*);
    ~Pentahedral();

    Element *clone() override;

    void renum(int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const;
    FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg=1) const;
    double getMass(const CoordSet& cs) const override;

    void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs);
    void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs=0);

    void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                     int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                      int surface=0, double *ndTemps=0);

    void markDofs(DofSetArray &) const override;
    int* dofs(DofSetArray &, int *p) const override;
     int numDofs() const override;

    int numNodes() const override;
    int* nodes(int * = 0) const override;

    int getTopNumber() override;

    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() const override { return 5; }
    int getDecFace(int iFace, int *fn);

    int getFace(int iFace, int *fn) { return getDecFace(iFace, fn); }

    void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
      { cCoefs = coefs; cFrame = frame; }

   double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
      { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                "              for Pentahedral el.\n"); return (double *) 0;
      }
    void getCFrame(CoordSet &cs, double cFrame[3][3]) const;

    void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                          int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

    void setMaterial(NLMaterial *) override;
    int numStates();
    void initStates(double *st);
    Corotator *getCorotator(CoordSet &cs, double *kel, int=2, int=2);
};

#endif
