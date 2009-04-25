#ifndef _NODAL_CONTACT_H_
#define _NODAL_CONTACT_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Element.d/Element.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/MyComplex.h>

//#include <Utils.d/list.h>
#include<list>
#include<map>
#include<utility>
using namespace std;

class PlanCommSetup {
public:
  PlanCommSetup() { initValue(); initTotalValue(); } //HB: more clean ..
  void setGlobStrNode(int GlobStrNode);
  void setGlobEndNode(int GlobEndNode);
  int getGlobStrNode();
  int getGlobEndNode();
  list<int> &getSubdList();
  void setCommunicate( bool Communicate);
  bool getCommunicate();
  void initValue();
  void addToValue( int component, double addValue);
  double *getValue();
  void initTotalValue();
  void addToTotalValue( double *addValue);
  double getTotalValue( int component);
  void initMultiplicity();
  void addToMultiplicity( int addValue);
  int getMultiplicity();
  void showData();
private:
  int globStrNode;
  int globEndNode;
  list<int> subdList;
  int lowestNormalDir;
  bool communicate;
  double value[3];
  double totalValue[3];
  int multiplicity; // number of Lagrange multiplier values that connect
                    // to this global contact
};

// (CKT) Karim TRAORE. May 2001
// Class dealing with unilateral or bilateral nodal contact

// typedef map<pair<int,int>, int> CtcMapType;
typedef map< pair<int,int>, list<int> > CtcMapType;
typedef CtcMapType::value_type CtcValuePair;

//--------------------------------------
class NodalContact {
  static double zeroForceTol;
  static bool haveDefaultFriction; // true if default friction coefficient has been given as input
  static bool haveFriction; // distinguish between problems having friction or no friction
  static int globalMode; // 0 => stitch OR 1 => contact
  int subNumber;
  PlanCommSetup *planCommSetup;
  int sPatSize;
  CtcMapType globCtcIdentMap;
  //list<int> **globCtcIdentLists;
public:
  void initPlanCommSetupPointer();
  PlanCommSetup &getPlanCommSetup(int pairNb);
  void allocatePlanCommSetup();
  void setPlanCommSetupGlobStrNodes(int *glNums, int numnodes);
  int &getSPatSize();
  void getGlobCtcIdent(int globStrNode, int globEndNode, list<int> &globCtcIdentList, int numnodes);
  // void makeGlobCtcIdent(int numnodes);
private:
  void findGlobCtcIdent(int GlobStrNode, int GlobEndNode,
                        list<int> &globCtcIdentList, int numnodes);
  int pairsNb;
  int pairsNbTotal;  // number of contact pairs for all subdomains
  int currentSize;
  int (*pairsList)[2]; // Gives (for a given local ctc numbering indice) the local sub nb of the couples of contact nodes
  // GR: new variable introduced to remove non-reliable comparison of double precision numbers in residual to zero
  bool *nullifyResidual; // used to record where the primal planning phase sets the residual to zero
  bool *computeGap; // used to indicate whether or not to compute the initial gap
  int *mode; // mode[iPair] = 0 => stitch OR mode[iPair] = 1 => contact
  double (*pair2Normal)[3];  //HB: store the current normal used in the node-to-node contact inequality
                             //    (set in setContactSurface or setContactNormal (see Salinas))
  // GR: (Tangent1,Tangent2,Normal) form in this order an orthonormal local Cartesian basis at contact node
  double (*pair2Tangent1)[3];
  double (*pair2Tangent2)[3];
  double (*pair2Tc1)[3]; // direction of frictional force
  double (*pair2Tc2)[3]; // normal to direction of frictional force in tangent plane
  bool (*pair2TgExists)[2];
  bool *pair2Tc1InActiveSet;
  bool *pair2Tc1NullifyPrimal;
  bool *pair2Tc2NullifyPrimal;
  bool *pair2Rconverged;
  int (*pair2SlidingStatus)[3];
  double *pair2SlipValue;
  static double commonCoulombFrictionCoef; // common value for the Coulomb friction coefficient
  double *pair2CoulombFrictionCoef; // storage for the Coulomb friction coefficient
  unsigned int *numConnectedNodes;
  int iEltStart ; // first virtual Elt nb (in domain)
  int iNodeStart; // first virtual Node nb (in domain)
  int (*pair2ccDofs)[3]; // pair2ccDofs[pair] contains the cc_dsa nb if not a corner. -1-(uc dof nb) otherwise. 
  int (*pair2cDofs)[3];  // pair2cDofs[pair] contains the c_dsa nb. -1 else
  int *virtualNode2Pair ; // Pair nb if virtual contact. -1 else
  double *ctcLMult; //HB: store the current value of the (normal) contact Lagrange multiplier
                    //    reset to 0.0 at every call to setContactSurface or setContactGap (for Salinas)
  double *gap0; //HB: store the (initial) rhs gap set in setContactSurface or setContactGap (for Salinas)
  double *gap;  //HB: store the current rhs gap used in the node-to-node contact inequality
  double *gapTg1;
  double *gapTg2;
  double *rotMat; 
  double ctcPrimalError; // the square of the ctc residual in the local sub
  double ctcError0; // initial displacement for lambda=0

  double epsilonPrimal;
  double epsilonDual;
  int *invBoundMap;
  int *invInternalMap;

  ResizeArray<double> *inactiveSetResidual; // store the negative part of the non active set residual 
                      // in order to update the residual from its old value (before the <>+ operation)
  int inactiveSetSize;
  ResizeArray<double> *inactiveSetFrictionResidual; // store the amount of slip that has been
                                                     // nullified in the jump (residual) in order 
                                                     // to update the residual from its old value
  int inactiveSetFrictionSize;

  int localNumNodes; 
  int globalNumNodes; // not including virtual nodes 
  bool *isCtcNode; // real or virtual node
  double *scaling;
  double fetiCtcTol;
  double *deltaFctc;
  double *deltaFctcTg1;
  double *deltaFctcTg2;
  ResizeArray<int> *dof2Pair; // boundary dof -> pairNb
  bool mergingNodes;
  double (*pair2rotMat)[9];
  bool flagCtc;
  int virtualElt(int i,int j);
  template <class Type, int d > int resize(Type (*&tab)[d], int dim, int checkDim);
  template <class Type> int resize(Type *&tab, int dim, int checkDim);

 public:
  int buildRotationMatrices();
  static const int virtualEltTypeNb=64;

  NodalContact();
  // GR: record subdomain number in NodalContact object
  NodalContact(NodalContact *, int SubNumber);
  ~NodalContact();
  void cleanAll(); //HB: delete all data
  void initializePointers();

  bool getFlagCtc();
  int addCtcPair(int n1, int n2, double nx, double ny, double nz,
                 double normalGap = 0.0, bool normalGapPresent = false,
                 double fricCoef = 0.0, bool fricCoefPresent = false,
                 int Mode = 1, bool modePresent = false);
  void setFlagCtc(bool flag);
  int getPairsNb();
  int getPairsNbTotal();
  double *getPair2Normal(int i);
  double *getPair2Tangent1(int i);
  double *getPair2Tangent2(int i);
  double *getPair2Tc1(int i);
  double *getPair2Tc2(int i);
  void setPair2Tc1InActiveSet(int i, bool active);
  void setPair2Tc1NullifyPrimal(int i, bool nullify);
  void setPair2Tc2NullifyPrimal(int i, bool nullify);
  void setPair2Rconverged(int i, bool converged);
  void setSlidingStatus(int i, int component, int setValue);
  int getSlidingStatus(int i, int component);
  void setSlipValue(int i, double slipValue);
  double getSlipValue(int i);
  bool isTc1InActiveSet(int i);
  bool isRconverged(int i);
  bool isTc1NullifiedInPrimal(int i);
  bool isTc2NullifiedInPrimal(int i);
  void constructTc2(int i); // use of this function assumes that Tc1 is known
  void reference2FrictionBasis(int i, double ug, double vg, double &uc, double &vc);
  void friction2ReferenceBasis(int i, double uc, double vc, double &ug, double &vg);
  bool * getPair2TgExists(int i);
  unsigned int &getNumConnectedNodes(int pairNb); // returns number of nodes 
                                                  // in this subdomain that
                                                  // share this pairNb number
  int getDirChangeRange();
  int getMode(int pairNb);
  int setPairsNb(int n, bool clean=false); //HB: flag to force to delete all arrays previously allocated
  int setPairsNbTotal(int n);
  int addVirtualEltInSub(Connectivity **stn, int numSub);
  int removeVirtualEltInSub(Connectivity **stn, int numSub);
  bool IsVirtualCtcNode(int NodeNb);
  bool IsCtcNode(int inode);
  bool IsCtcYes(int pairNb); // DEBUG
  int makeIsCtcNode(int numnodes);
  int makeIsCtcNodeSalinas(int numnodes);
  int getVirtualNodeNb(int NodeNb);
  int renum(int *table, int glmax, int numnodes, int *glNums);
  int *getPair2ccDofs(int nbPairs);
  int *getPair2cDofs(int nbPairs);
  int setPair2Dofs(int i, int *dofcc, int *dofc);
  void changeSign(int pair);
  void createNonSharedCtcNumbers(int sharedDim, int numNodes);
  int sign(int pair); // subdomain level
  double dsign(int pair);
  int computeInitialGap(CoordSet &nodes);
  double setNormalGap(int pairNb, double ngap);  //HB: gap = gap0 <- ngap
  double updateNormalGap(int pairNb, double Dg, bool first); //HB: gap <- gap0 - Dg
  void setNormal(int pairNb, double* nrm);       //HB
  double& viewContactMuliplier(int pairNb) { return(ctcLMult[pairNb]); }; //HB
  void zeroCtcLagrangeMult(); //HB
  double &Gap(int pairNb, int pairTg = 0);
  bool IsInActiveSet(int pair);
  int changeStatus(int pair);
  int computeCtcPrimalError(int *allBoundDofs, int interfSize, double *u, int *boundDofFlag); // PJSA
  int computeCtcPrimalError(int *allBoundDofs, int interfSize, DComplex *u, int *boundDofFlag) { return 0; };
  void setCtcPrimalError(double error);
  double getCtcPrimalError();
  double getEpsilonPrimal();
  double getEpsilonDual();
  double getCtcError0();
  void setCtcError0(double rhs);
  void setInvBoundMap(int len, int *tab);
  void setInvInternalMap(int len, int *tab);
  int *getInvBoundMap();
  int *getInvInternalMap();
  void initInactiveSetResidual();
  double saveAndCheckInactiveSetResidual(double value);
  double getInactiveSetResidual();
  void initInactiveSetFrictionResidual();
  void saveInactiveSetFrictionResidual( double Tg1Value, double Tg2Value);
  double getInactiveSetFrictionResidual();
  double getPair2Scaling(int pairNb);
  double getFetiCtcTol();
  void setDeltaFctc(int pairNb, int pairTg, double val);
  double getDeltaFctc(int pairNb, int pairTg);
  void setDof2Pair(int nodeNb,int dofNb,int nnn); 
  int getDof2Pair(int dofNb);
  void setNodeMerging(int pair);
  bool nodeMerging();
  void nodeMerging(int num, int etype, int nnodes, int *n);
  int rotateTranspose(double *u, double *tRu, int i);
  int rotate(double *u, double *Ru, int i);

  void resetNullifyResidual();
  void setNullifyResidual(int pair, bool setValue);  // replaces activateNullifyResidual(int pair)
  bool isNullifiedInResidual(int pair);

  bool hasFriction();
  void setCommonCoulombFricCoef(double CommonCoulombFricCoef);
  void setGlobalMode(int glbMode);
  double &getPair2CoulombFricCoef(int pairNb);
  double getZeroForceTol();
  void setUnitTangents();

  int *operator[] (int i);

  // preprocessing
  int setRotMat(double d0, double d1, double d2,
         	double d3, double d4, double d5,
		double d6, double d7, double d8,
		double d9, double d10, double d11);
  int rotate(double **M); 
  int rotate(double *Mx,double *My,double *Mz); 
  bool rotationDefined();

  // initialization for salinas
  void initSandia(int iNodeStart, int numNodes, int *glNums); 
  void initSandia();
  void computeScaling(int *_nL, int *_nR, int *nodeToSubCount);

  double *getNormalOrTangent(int pairNb, int pairTg);

  void print(int* subNum=0);
};

#endif
