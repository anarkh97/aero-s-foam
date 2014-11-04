#ifndef _SUPERELEMENT_H_
#define _SUPERELEMENT_H_

#include <Element.d/Element.h>
#include <map>
class SuperCorotator;

class SuperElement : public Element 
{
  private:
    Elemset *eset;
    DofSetArray *dsa;
    SuperCorotator *superCorotator;
    int nInternalNodes;
    double **sub_extf;

  protected:
    CoordSet *css;
    Element **subElems;    
    int nSubElems;       
    int **subElemDofs;
    int **subElemNodes;
    int nnodes; // not including internal nodes
    int ndofs;
    int *nn; // all the node numbers
    bool localFlag;

    FullSquareMatrix stiffness(CoordSet& cs, double *k, int flg=1);
  public:
    SuperElement(bool = false);
    virtual ~SuperElement();

    double * getPreviouslyComputedSubExternalForce(int i) { return (sub_extf) ? sub_extf[i] : 0; }

    int getNumSubElems() { return nSubElems; }
    int getSubElemNumDofs(int i) { return subElems[i]->numDofs(); }
    int getSubElemNumNodes(int i) { return subElems[i]->numNodes(); }
    int* getSubElemDofs(int i) { return subElemDofs[i]; }
    int* getSubElemNodes(int i) { return subElemNodes[i]; }

    void setPressure(PressureBCond *);
    PressureBCond* getPressure();

    void renum(int *table);
    void renum(EleRenumMap&);
    void setGlNum(int gn, int sn = 0);

    void setProp(StructProp *p, bool _myProp = false); 
    void setPreLoad(std::vector<double> &load);
    std::vector<double> getPreLoad();
    void setFrame(EFrame *frame);
    void buildFrame(CoordSet &cs);
    void setOffset(double *o);
    void setCompositeData(int _type, int nlays, double *lData,
                          double *coefs, double *frame);
    double * setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta);
    void setMaterial(NLMaterial *);

    FullSquareMatrix massMatrix(CoordSet& cs, double *m, int cmflg=1);
    void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg=1);
    void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);

    double getMass(CoordSet&);
    double getMassThicknessSensitivity(CoordSet&);
    double weight(CoordSet&, double *);
    double getWeightThicknessSensitivity(CoordSet&, double *);
    void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration);
    void getGravityForce(CoordSet &cs, double *gravity, Vector &force,
                         int gravflg, GeomState *gs=0);
    void getGravityForceThicknessSensitivity(CoordSet &cs, double *gravity, Vector &forceSen,
                                             int gravflg, GeomState *gs=0);
    void getGravityForceNodalCoordinateSensitivity(CoordSet &cs, double *gravity, GenFullM<double> &,
                                                   int gravflg, GeomState *gs=0);
    void getThermalForce(CoordSet &cs, Vector &ndT, Vector &force,
                         int glflag, GeomState *gs=0);
    void getIntrnForce(Vector &elForce, CoordSet &cs,
                       double *elDisp, int Index, double *ndTemps);
    void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                     Vector &elDisp, int strInd, int surface=0,
                     double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);
    void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                         double *, int avgnum, double ylayer, double zlayer);
    void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                            double *ndTemps, int avgnum, double ylayer, double zlayer);
    void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                               double * = 0, int avgnum=1, double ylayer=0, double zlayer=0);
    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                      Vector &elDisp, int strInd, int surface=0,
                      double *ndTemps=0);
    void computeHeatFluxes(Vector &heatflux, CoordSet &cs, Vector &elTemp, int hflInd);
    void trussHeatFluxes(double &trussflux, CoordSet &cs, Vector &elTemp, int hflInd);
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *, double *res, GeomState *gs=0);
    void computeTemp(CoordSet &cs, State &state, double[2], double *res);
    void getFlFlux(double[2], double *flF, double *res);

    void markDofs(DofSetArray &dsa);
    int* dofs(DofSetArray &dsa, int *p=0);
    int numDofs();
    int numNodes();
    int* nodes(int *p=0);

    Corotator *getCorotator(CoordSet &, double *, int = 2, int = 2);
    void computePressureForce(CoordSet& cs, Vector& elPressureForce,
                              GeomState *gs=0, int cflg = 0, double t = 0);

    double* getMidPoint(CoordSet &cs);
    double* getCompositeData(int nl);
    double* getCompositeFrame();
    int getCompositeLayer();
    int dim();
    void addFaces(PolygonSet *pset);
    int numInternalNodes();
    void setInternalNodes(int *in);
    bool isSafe();
    bool isRotMidSideNode(int iNode);
    bool isMpcElement();
    //bool isRigidMpcElement(const DofSet & = DofSet::nullDofset, bool forAllNodes=false);
    bool isConstraintElement();

    int getMassType();
    int getNumMPCs();
    LMPCons** getMPCs();
    void makeAllDOFs();

    int numStates();
    void setStateOffset(int);
    void initStates(double *);
};

#endif
