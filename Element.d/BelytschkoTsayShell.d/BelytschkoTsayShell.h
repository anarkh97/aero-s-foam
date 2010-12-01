#ifndef _BELYTSCHKOTSAYSHELL_H_
#define _BELYTSCHKOTSAYSHELL_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

class GeomState;
class MultiFront;
class NLMaterial;
class ExpMat;

class BelytschkoTsayShell : virtual public Element, public Corotator
{
  public:
    static double t1, t2, t3, t4, t5, t6, t7;
  protected:
    // TODO most of this should belong to element property (and therefore be shared to save memory)
    int nn[4];
    int optele; // element type option (3 for bt shell)
    int optmhd; // numerical method option (0 for conventional fem)
    int optdmg; // damage model type (0 for no damage, 1 for lematire damage model, 2 for linear softening with scaling)
    int opthgc; // hourglass control (1 for perturbation type hourglass control)
    int optcri[2]; // crack criterion
    int opttrc; // bc option (0 for pressure, 1 for traction, -1 for neither)
    int optdmp; // damping (0/1 for damping off/on)
    int optcor[2]; // warping and/or shear correction on/off
    double prmhgc[10]; // hourglass control parameters
    double prmdmp[10]; // damping control parameters
    int ngqpt[3]; // ngqpt[0] = gq rule for regular element
                  // ngqpt[1] = gq rule for enriched element
                  // ngqpt[3] = gq rule for through thickness
    int ngqpt4; // gq rule for bc or cohesive force integration
    int nndof;
    int ndime;
    int nnode;

    int mgaus[3];
    int mgqpt[2];
    double *gqpoin3;
    double *gqweigt3;
    double *evar1; // effective strain and damage
    double *evar2; // effective stress
    double *evoit1; // voight form of hourglass control stress
    double *evoit2; // voight form of local cauchy stress
    double *evoit3; // strain (local)
    ExpMat *expmat;
    MFTTData *mftt;

  public:
    BelytschkoTsayShell(int*);
    ~BelytschkoTsayShell();

    //void setProp(StructProp *p, bool _myProp = false);
    void setMaterial(NLMaterial *);
    void setPressure(double, MFTTData* = 0);
    double getPressure();
    Element *clone();

    void renum(int *);

    FullSquareMatrix stiffness(CoordSet&, double* d, int flg = 1);
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
                         GeomState *gs);
    void getVonMises(Vector& stress, Vector& weight, CoordSet& cs, 
                     Vector& elDisp,  int strInd, int surface = 0,
                     double *ndTemps = 0, double ylayer = 0.0,
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
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double delt, double);

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

    void writeHistory(int fn);
    void readHistory(int fn);

};
#endif

