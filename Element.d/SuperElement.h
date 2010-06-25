#ifndef _SUPERELEMENT_H_
#define _SUPERELEMENT_H_

#include <Element.d/Element.h>
#include <map>
class SuperCorotator;

class SuperElement : public Element 
{
  protected:
        // all of the following member variables need to be set in the constructor
        Elemset *eset;
        DofSetArray *dsa;
        Element **subElems;    
        int nSubElems;       
        SuperCorotator *superCorotator;
        int **subElemDofs;
        int **subElemNodes;
        int nnodes;
        int ndofs;
        int *nn; // all the node numbers

        FullSquareMatrix stiffness(CoordSet& cs, double *k, int flg=1);
  public:
        SuperElement() { nn = 0; superCorotator = 0; };
        virtual ~SuperElement();

        int getNumSubElems() { return nSubElems; }
        int getSubElemNumDofs(int i) { return subElems[i]->numDofs(); }
        int getSubElemNumNodes(int i) { return subElems[i]->numNodes(); }
        int* getSubElemDofs(int i) { return subElemDofs[i]; }
        int* getSubElemNodes(int i) { return subElemNodes[i]; }

        void setPressure(double pres);
        double getPressure() { return subElems[0]->getPressure(); }

        void renum(int *table);

        void setProp(StructProp *p, bool _myProp = false); 
        void setPreLoad(double load, int &flg);
        void setFrame(EFrame *frame);
        void buildFrame(CoordSet &cs);
        void setOffset(double *o);
        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame);
        double * setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta);

        FullSquareMatrix massMatrix(CoordSet& cs, double *m, int cmflg=1);

        double getMass(CoordSet&);
        void getGravityForce(CoordSet &cs, double *gravity, Vector &force,
                             int gravflg, GeomState *gs=0);
        void getThermalForce(CoordSet &cs, Vector &ndT, Vector &force,
                             int glflag, GeomState *gs=0);
        void getIntrnForce(Vector &elForce, CoordSet &cs,
                           double *elDisp, int Index, double *ndTemps);
        void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface=0,
                         double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);
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
                                  GeomState *gs=0, int cflg = 0);

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

        void initialize(int, int*);
};

#endif
