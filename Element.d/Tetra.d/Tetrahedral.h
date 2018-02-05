#ifndef _TETRAHEDRAL_H_
#define _TETRAHEDRAL_H_

#include <Element.d/Element.h>
#include <Element.d/Tetra.d/TetraElementTemplate.hpp>

class Tetrahedral: public Element,
                   public TetraElementTemplate<double>
{
    int nn[4];
    double *cCoefs;
    double *cFrame;
    NLMaterial *mat;
    void computeDjDx(double x[4], double y[4], double z[4], double J, double djdx[12]);

  public:
    Tetrahedral(int*);
    ~Tetrahedral();

    Element *clone() override;

    void renum(int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const;
    void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
    FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
    void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration);
    double getMass(const CoordSet& cs) const;
    void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega);

    void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs);
    void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                   GenFullM<double> &dGfdx, int gravflg, GeomState *geomState);
    void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs=0);

    void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                     int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight,
                                               CoordSet &cs, Vector &elDisp, int strInd,
                                               int surface, double* ndTemps=0,
                                               int avgnum=1, double ylayer=0, double zlayer=0);

    void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *,
                                            CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                            double *ndTemps, int avgnum, double ylayer, double zlayer);

    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                      int surface=0, double *ndTemps=0);

    void markDofs(DofSetArray &) override;
    int* dofs(DofSetArray &, int *p) override;
     int numDofs() const override;

    int numNodes() const override;
    int* nodes(int * = 0) const override;

    int getTopNumber() override;

    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() const override { return 4; }
    int getDecFace(int iFace, int *fn);

    int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }

    void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
      { cCoefs = coefs; cFrame = frame; }

   double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
      { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                "              for Tetrahedral el.\n"); return (double *) 0;
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
