#ifndef _BELYTSCHKOTSAYSHELL_H_
#define _BELYTSCHKOTSAYSHELL_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

class GeomState;
class MultiFront;

class BelytschkoTsayShell : virtual public Element, public Corotator
{
  protected:
    int nn[4];
    int optele; // element type option (3 for bt shell)
    int optmhd; // numerical method option (0 for conventional fem)
    int optctv; // constiutive law (1 for hypoelastic, 3 for elasto viscoplastic, 5 for j2 explicit)
    int optdmg; // damage model type (0 for no damage, 1 for lematire damage model, 2 for linear softening with scaling)
    int opthgc; // hourglass control (1 for perturbation type hourglass control)
    int optcri[10]; // crack criterion
    double prmdmg[10]; // damage control parameters
    double prmhgc[10]; // hourglass control parameters
    int ngqpt[3]; // ngqpt[0] = gq rule for regular element
                  // ngqpt[1] = gq rule for enriched element
                  // ngqpt[3] = gq rule for through thickness
    int nndof;
    int ndime;
    int nnode;

  public:
    BelytschkoTsayShell(int*);

    Element *clone();

    void renum(int *);

    FullSquareMatrix stiffness(CoordSet&, double* d, int flg = 1);
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
                         GeomState *gs);
    void getVonMises(Vector& stress, Vector& weight, CoordSet& cs, 
                     Vector& elDisp, int strInd, int surface = 0,
                     double *ndTemps = 0,double ylayer = 0.0,
                     double zlayer = 0.0, int avgnum = 0);
    void getAllStress(FullM& stress, Vector& weight, CoordSet& cs,
                      Vector& elDisp, int strInd, int surface = 0,
                      double* ndTemps = 0);

    void markDofs(DofSetArray&);
    int* dofs(DofSetArray&, int* p = 0);
    int numDofs();

    int numNodes();
    int* nodes(int* = 0);
    Corotator *getCorotator(CoordSet&, double*, int , int);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*);

    void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
                    GeomState*);
    void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
                   GeomState* = 0);

    int getTopNumber();
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState* gs, int cflg);
              
    void getThermalForce(CoordSet& cs, Vector& ndTemps, Vector &elThermalForce, 
                     int glfag, GeomState* gs = 0);
                                        
    bool isShell() { return true; }

    int getMassType() { return 0; } // lumped only

    bool hasRot() { return true; }

    PrioInfo examine(int sub, MultiFront* mf);
};
#endif

