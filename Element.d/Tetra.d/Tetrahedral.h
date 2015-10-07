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

    Element *clone();

    void renum(int *);
    void renum(EleRenumMap&);

    FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
    void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
    void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration);
    double getMass(CoordSet& cs);
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

    void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                            double *ndTemps, int avgnum, double ylayer, double zlayer);

    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                      int surface=0, double *ndTemps=0);

    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int numDofs();

    int numNodes();
    int* nodes(int * = 0);

    int getTopNumber();

    PrioInfo examine(int sub, MultiFront *);
    int nDecFaces() { return 4; }
    int getDecFace(int iFace, int *fn);

    int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }

    void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame)
      { cCoefs = coefs; cFrame = frame; }

    double* setCompositeData2(int, int, double*, double*, CoordSet&, double)
      { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                "              for Tetrahedral el.\n"); return (double *) 0;
      }

    void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                          int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

    void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

    void setMaterial(NLMaterial *);
    int numStates();
    Corotator *getCorotator(CoordSet &cs, double *kel, int=2, int=2);
};

#endif
