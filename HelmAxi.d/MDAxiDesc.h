#ifndef _MDAXIDESC_H_
#define _MDAXIDESC_H_

template<class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
class Connectivity;
class DistrComplexVector;
class Domain;
class FetiHAxiSolver;
class FourierHelmBCs;
class MDAxiData;
class MPCData;
class ScatterData;
class SommerElement;
class StaticTimers;


class MDAxiDesc {

     Domain *domain;

     MDAxiData **localData;
     int numSub;

     Connectivity *subToElem, *subToNode, *nodeToSub, *subToSub;

     FourierHelmBCs *glBCs;
     ScatterData *glScatter;

     FetiHAxiSolver *solver;

     DistrComplexVector *sol;
     DistrComplexVector *rhs;

     MPCData *glMPCs;
     Connectivity *mpcToNode, *mpcToSub;

     StaticTimers *times;

     public :

     MDAxiDesc(Domain *d, FourierHelmBCs *fbcs, MPCData *mpcs, ScatterData *);

     void preProcess();
     void getDecomp(FILE *);
     void makeSubD();
     void buildMpcToNode();
     void constructSubDomains(int iSub, MDAxiData **ld, 
                         Domain *d, Connectivity *cn, Connectivity *sToN);
     void getDomainConnect();
     void distributeBCHelm(int);
     void distributeMPCs();
     void extractSubDomainMPCs(int);
     void buildInterfacePolygons();
     void finishInterfacePolygons();

     void makeInterface();
     void prepareCoarseData(int, DofSet ***, int *, int **, int **, int **,
                            DComplex ***);
     void makeCoarseData(int, DofSet ***, int *, int **, int **, int **,
                            DComplex ***);
     void deleteCoarseInfo(int, int **, int **, int **, DComplex ***);
     void prepareMPCSet(int );
 
     void getInterfSigns();
     void solve();
     void buildRhsMPC(ComplexVector &);
     void postProcessing(DistrComplexVector &);

     void buildFFP(DComplex **, DistrComplexVector &);

     void termDUDN(int nT, DComplex *FFP, DComplex* *dudn, 
                   double (*vectorDir)[3], int *pointNumb);
     void ffpAxiNeum(int, DComplex *, DComplex**, double (*)[3]);
     void ffpAxiDir(int, DComplex *, DComplex* *, double (*)[3]);

     void assembleDUDN(int, DistrComplexVector *, DComplex **);
     void localDUDN(int, DComplex **, DComplex **, int *, int *);

     void verifMPCScat();

};

#endif
