#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <Utils.d/dofset.h>
#include <Contact.d/NodalContact.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/Domain.h>
#include <Utils.d/b2Util.h>
#include <Utils.d/b2Rotation.h>

extern int verboseFlag;

// (CKT) Karim TRAORE. May 2001
// Class dealing with unilateral or bilateral nodal contact

double NodalContact::zeroForceTol = 1.0e-12;
bool NodalContact::haveDefaultFriction = false;
bool NodalContact::haveFriction = false;
int NodalContact::globalMode = 1;
double NodalContact::commonCoulombFrictionCoef = 0.0;

NodalContact::NodalContact()  
{
  initializePointers();
  currentSize = 0;
  pairsNb = 0;
  pairsNbTotal = 0;
  sPatSize = 0; 
  subNumber = -1;
  initPlanCommSetupPointer();
  iEltStart = -1;
  iNodeStart = -1;
  ctcPrimalError = 0;
  ctcError0 = 1;
  epsilonPrimal = 1E-4;
  //  GR: epsilonDual reset to 1 since the limiting value for planning
  //  actions during dual planning is now implemented such that it
  //  will decrease proportional to the precondition residual
  epsilonDual = 1.0;
  fetiCtcTol = 1E-6;
  inactiveSetSize = 0;
  inactiveSetFrictionSize = 0;
  mergingNodes = false;
  pair2rotMat = 0;
  flagCtc = false;
  globalNumNodes = 0;
};

//------------------------------------------------------
NodalContact::NodalContact(NodalContact *ctcDom, int SubNumber) 
{
  initializePointers();
  currentSize = ctcDom->currentSize;
  sPatSize = 0;  
  pairsNb = ctcDom->pairsNb;
  pairsList = ctcDom->pairsList;
  computeGap = ctcDom->computeGap;
  mode = ctcDom->mode;
  subNumber = SubNumber;
  pair2CoulombFrictionCoef = ctcDom->pair2CoulombFrictionCoef;
  pair2Tangent1 = ctcDom->pair2Tangent1;
  pair2Tangent2 = ctcDom->pair2Tangent2;
  pair2Tc1 = ctcDom->pair2Tc1;
  pair2Tc2 = ctcDom->pair2Tc2;
  pair2TgExists = ctcDom->pair2TgExists;
  pair2Tc1InActiveSet = ctcDom->pair2Tc1InActiveSet;
  pair2Tc1NullifyPrimal = ctcDom->pair2Tc1NullifyPrimal;
  pair2Tc2NullifyPrimal = ctcDom->pair2Tc2NullifyPrimal;
  gapTg1 = ctcDom->gapTg1;
  gapTg2 = ctcDom->gapTg2;
  initPlanCommSetupPointer();
  pair2Normal=ctcDom->pair2Normal;  
  pair2ccDofs = ctcDom->pair2ccDofs;
  pair2cDofs = ctcDom->pair2cDofs;
  virtualNode2Pair = ctcDom->virtualNode2Pair;
  gap0 = ctcDom->gap;
  gap  = ctcDom->gap;
  iEltStart = ctcDom->iEltStart;
  iNodeStart = ctcDom->iNodeStart;
  isCtcNode = ctcDom->isCtcNode;
  ctcPrimalError = 0;
  ctcError0 = 1;
  epsilonPrimal = ctcDom->epsilonPrimal;
  epsilonDual = ctcDom->epsilonDual;
  inactiveSetSize = 0;
  inactiveSetFrictionSize = 0;
  scaling = ctcDom->scaling;
  mergingNodes = false;
  flagCtc = false;
  globalNumNodes = 0;
}
  
/***************************************************************************/
// initialization for salinas, PJSA 6-28-02
// since ctcDom is never built and renum is not called  

#define HB_CTC_DEBUG
void NodalContact::initSandia(int iNodeStart, int numnodes, int *glNums) 
{  
  int i;

  pair2ccDofs=new int[pairsNb][3];
  pair2cDofs=new int[pairsNb][3];
  deltaFctc = new double [pairsNb];

  // build virtualNode2Pair, PJSA: modified for MPCs 12-19-02
  virtualNode2Pair = new int [numnodes];
  for(i=0; i<numnodes; ++i) virtualNode2Pair[i] = -1;
#ifdef HB_CTC_DEBUG
  if(iNodeStart+pairsNb!=numnodes)
    fprintf(stderr," *** WARNING: problem in NodalContact::initSandia: iNodeStart = %d, numnodes = %d, pairsNb = %d, iNodeStart+pairsNb = %d\n");
#endif
  for (i=iNodeStart; i<(iNodeStart+pairsNb); ++i) 
    virtualNode2Pair[i] = i - iNodeStart;  // pair number

  if(!numConnectedNodes) { //HB 06-23-05
    numConnectedNodes = new unsigned int[pairsNb];
    for(int iPair = 0; iPair < pairsNb; iPair++)
      numConnectedNodes[iPair] = 0;
  }

  if(!pair2SlidingStatus) { //HB 06-23-05
    pair2SlidingStatus = new int[pairsNb][3];
    int kPair;
    for(int iPair = 0; iPair < pairsNb; iPair++) 
      this->setSlidingStatus(iPair, 1, 0);
  }
  
  if(!ctcLMult) {
    ctcLMult = new double[pairsNb];
    zeroCtcLagrangeMult();
  }
}

void NodalContact::zeroCtcLagrangeMult()
{
  if(!ctcLMult) ctcLMult = new double[pairsNb];
  for(int iPair = 0; iPair < pairsNb; iPair++)
      ctcLMult[iPair] = 0.0; 
}

void NodalContact::initSandia()
{
  // PJSA 11-3-04 version for new CU_Feti wrapper
  isCtcNode = 0; // shouldn't be used anymore

  pair2ccDofs=new int[pairsNb][3];
  pair2cDofs=new int[pairsNb][3];
  deltaFctc = new double [pairsNb];

  // virtual node number now is the same as the pair number
  // and should eventually be eliminated
  virtualNode2Pair = new int [pairsNb];
  for(int i=0; i<pairsNb; ++i) virtualNode2Pair[i] = i;

  buildRotationMatrices();
  initPlanCommSetupPointer();
  sPatSize = pairsNb;
                                                                                                                             
  nullifyResidual = new bool[pairsNb];
  resetNullifyResidual();

  if(!numConnectedNodes) {
    numConnectedNodes = new unsigned int[pairsNb];
    for(int i = 0; i < pairsNb; i++) numConnectedNodes[i] = 0;
  }

  if(!pair2SlidingStatus) {
    pair2SlidingStatus = new int[pairsNb][3];
    for(int i = 0; i < pairsNb; i++) setSlidingStatus(i, 1, 0);
  }

  if(!ctcLMult) {
    ctcLMult = new double[pairsNb];
    zeroCtcLagrangeMult();
  }
}

void NodalContact::computeScaling(int *_nL, int *_nR, int *nodeToSubCount)
{
  int i;
  for(i=0; i<pairsNb; i++) {
    int nL = _nL[i];
    int nR = _nR[i];
    scaling[i] = 1.0E+00 / (nodeToSubCount[nL] + nodeToSubCount[nR]);
  }
//  std::cerr << "scaling = "; for(i=0; i<pairsNb; i++) std::cerr << scaling[i] << " "; std::cerr << std::endl;
}

/****************************************************************************/

//------------------------------------------------------
NodalContact::~NodalContact() 
{
/*  if(planCommSetup) delete [] planCommSetup;

  if(haveFriction) {
    if(pair2Tangent1) delete [] pair2Tangent1;
    if(pair2Tangent2) delete [] pair2Tangent2;
    if(pair2Tc1) delete [] pair2Tc1;
    if(pair2Tc2) delete [] pair2Tc2;
    if(pair2TgExists) delete [] pair2TgExists;
    if(pair2Tc1InActiveSet) delete [] pair2Tc1InActiveSet;
    if(pair2Tc1NullifyPrimal) delete [] pair2Tc1NullifyPrimal;
    if(pair2Tc2NullifyPrimal) delete [] pair2Tc2NullifyPrimal;
    if(pair2Rconverged) delete [] pair2Rconverged;
    if(pair2SlipValue) delete [] pair2SlipValue;
    if(gapTg1) delete [] gapTg1;
    if(gapTg2) delete [] gapTg2;
  }

  if(pairsList) delete [] pairsList;
  if(pair2Normal) delete [] pair2Normal;
  if(gap) delete [] gap;
  if(scaling) delete [] scaling;
  if(pair2CoulombFrictionCoef) delete [] pair2CoulombFrictionCoef;
  if(computeGap) delete [] computeGap;
  if(mode) delete [] mode;
  if(pair2SlidingStatus) delete [] pair2SlidingStatus; 
  if(isCtcNode) delete [] isCtcNode;
  if(virtualNode2Pair) delete [] virtualNode2Pair;

  if(nullifyResidual) delete [] nullifyResidual; 
  if(numConnectedNodes) delete [] numConnectedNodes;
  if(pair2ccDofs) delete [] pair2ccDofs;
  if(pair2cDofs) delete [] pair2cDofs;
  if(rotMat) delete [] rotMat;
  if(invBoundMap) delete [] invBoundMap;
  if(invInternalMap) delete [] invInternalMap;
  if(inactiveSetResidual) delete inactiveSetResidual;
  if(inactiveSetFrictionResidual) delete inactiveSetFrictionResidual;
  if(deltaFctc) delete [] deltaFctc;
  if(deltaFctcTg1) delete [] deltaFctcTg1;
  if(deltaFctcTg2) delete [] deltaFctcTg2;
  if(dof2Pair) delete dof2Pair;
  if(pair2rotMat) delete [] pair2rotMat;
*/
  cleanAll();
}

void
NodalContact::cleanAll()
{
  if(planCommSetup) delete [] planCommSetup;

  if(haveFriction) {
    if(pair2Tangent1) delete [] pair2Tangent1;
    if(pair2Tangent2) delete [] pair2Tangent2;
    if(pair2Tc1) delete [] pair2Tc1;
    if(pair2Tc2) delete [] pair2Tc2;
    if(pair2TgExists) delete [] pair2TgExists;
    if(pair2Tc1InActiveSet) delete [] pair2Tc1InActiveSet;
    if(pair2Tc1NullifyPrimal) delete [] pair2Tc1NullifyPrimal;
    if(pair2Tc2NullifyPrimal) delete [] pair2Tc2NullifyPrimal;
    if(pair2Rconverged) delete [] pair2Rconverged;
    if(pair2SlipValue) delete [] pair2SlipValue;
    if(gapTg1) delete [] gapTg1;
    if(gapTg2) delete [] gapTg2;
  }

  if(pairsList) delete [] pairsList;
  if(pair2Normal) delete [] pair2Normal;
  if(gap0) delete [] gap0;
  if(gap) delete [] gap;
  if(scaling) delete [] scaling;
  if(pair2CoulombFrictionCoef) delete [] pair2CoulombFrictionCoef;
  if(computeGap) delete [] computeGap;
  if(mode) delete [] mode;
  if(pair2SlidingStatus) delete [] pair2SlidingStatus;
  if(isCtcNode) delete [] isCtcNode;
  if(virtualNode2Pair) delete [] virtualNode2Pair;

  if(nullifyResidual) delete [] nullifyResidual;
  if(numConnectedNodes) delete [] numConnectedNodes;
  if(pair2ccDofs) delete [] pair2ccDofs;
  if(pair2cDofs) delete [] pair2cDofs;
  if(rotMat) delete [] rotMat;
  if(invBoundMap) delete [] invBoundMap;
  if(invInternalMap) delete [] invInternalMap;
  if(inactiveSetResidual) delete inactiveSetResidual;
  if(inactiveSetFrictionResidual) delete inactiveSetFrictionResidual;
  if(deltaFctc) delete [] deltaFctc;
  if(deltaFctcTg1) delete [] deltaFctcTg1;
  if(deltaFctcTg2) delete [] deltaFctcTg2;
  if(dof2Pair) delete dof2Pair;
  if(pair2rotMat) delete [] pair2rotMat;

  if(ctcLMult) delete [] ctcLMult;

  initializePointers();
}

void
NodalContact::initializePointers()
{
  planCommSetup = 0; pairsList = 0; nullifyResidual = 0; computeGap = 0; mode = 0;
  pair2Normal = 0; pair2Tangent1 = 0; pair2Tangent2 = 0; pair2Tc1 = 0; pair2Tc2 = 0;
  pair2TgExists = 0; pair2Tc1InActiveSet = 0; pair2Tc1NullifyPrimal = 0;
  pair2Tc2NullifyPrimal = 0; pair2Rconverged = 0; pair2SlidingStatus = 0;
  pair2SlipValue = 0; pair2CoulombFrictionCoef = 0; numConnectedNodes = 0;
  pair2ccDofs = 0; pair2cDofs = 0; virtualNode2Pair = 0; gap0 = 0; gap = 0; gapTg1 = 0;
  gapTg2 = 0; rotMat = 0; invBoundMap = 0; invInternalMap = 0; inactiveSetResidual = 0;
  inactiveSetFrictionResidual = 0; isCtcNode = 0; scaling = 0; deltaFctc = 0; 
  deltaFctcTg1 = 0; deltaFctcTg2 = 0; dof2Pair = 0; pair2rotMat = 0;
  ctcLMult = 0; 
}


//------------------------------------------------------
int
NodalContact::setPairsNb(int n, bool clean) 
{
  if(clean) cleanAll(); //HB: check for memory leaks when used with Salinas 
  if((pairsNb != 0) || (n < 1)) {
    if(verboseFlag) std::cout << "__ Invalid contact pairs number prescription" << std::endl;
    exit(-1);
  }
  
  pairsNb                  = n;
  pairsList                = new int[n][2];
  pair2CoulombFrictionCoef = new double[n];
  computeGap               = new bool[n];
  mode                     = new int[n];
  pair2Normal              = new double[n][3];
  gap0                     = new double[n];
  gap                      = new double[n];
  pair2ccDofs              = 0;
  pair2cDofs               = 0;
  scaling                  = new double[n];
  ctcLMult                 = new double[n];
  return 0;
}

// GR ------------------------------------------------------
void
NodalContact::setCommonCoulombFricCoef( double CommonCoulombFricCoef)
{
  commonCoulombFrictionCoef = CommonCoulombFricCoef;
  haveDefaultFriction = true;
  // globalMode = 1;
}

void
NodalContact::setGlobalMode( int glbMode)
{
  if(glbMode == 0) {
    globalMode = 0;
    // commonCoulombFrictionCoef = 1.0;
    // haveDefaultFriction = true;
  }
  else {
    globalMode = 1;
    // commonCoulombFrictionCoef = 0.0;
    // haveDefaultFriction = true;
  }
}

// GR --------------------------------------------------
void
NodalContact::resetNullifyResidual()
{
  if(!nullifyResidual) { nullifyResidual = new bool[pairsNb]; } //HB
  for(unsigned int i = 0; i < pairsNb; i++) nullifyResidual[i] = false;
}

//------------------------------------------------------
int 
NodalContact::makeIsCtcNode(int numnodes) 
{ 
  if(currentSize < pairsNb) {
    filePrint(stderr, "*** ERROR: Less contact pairs then expected ***\n");
    exit(-1); // can not define more than the described number of contact nodes
  }

  if(isCtcNode) delete isCtcNode;
  isCtcNode = new bool [numnodes];  // PA: bug fix
  globalNumNodes = numnodes;

  int i;
  for(i=0; i<numnodes; i++) 
    isCtcNode[i] = false;

  for(i=0; i<pairsNb; i++) {
    isCtcNode[pairsList[i][0]] = true; 
    isCtcNode[pairsList[i][1]] = true;
  }

  int count = 0;
  for(i=0; i<numnodes; i++) {
    if(isCtcNode[i]) count++;
  }

  return 0;
}

// GR:
int 
NodalContact::makeIsCtcNodeSalinas(int numnodes) 
{ 
  int i;
  isCtcNode = new bool [numnodes];
  
  for(i = 0; i < numnodes; i++) 
    isCtcNode[i] = false;

  for(i = 0; i < pairsNb; i++)
    isCtcNode[pairsList[i][0]] = true;

  return 0;
}

//------------------------------------------------------
int
NodalContact::addCtcPair(int n1, int n2, double nx, double ny, double nz,
                         double normalGap, bool normalGapPresent,
                         double fricCoef, bool fricCoefPresent,
                         int Mode, bool modePresent) 
{
  // static int currentSize = 0;  PJSA -> made this a class variable
  if((currentSize + 1) > pairsNb) {
    filePrint(stderr, "*** ERROR: More contact pairs then expected ***\n"); 
    exit(-1); // can not define more than the described number of contact nodes
  }

  if (rotationDefined()) rotate(&nx, &ny, &nz);

  pairsList[currentSize][0]=n1;
  pairsList[currentSize][1]=n2;
  // GR: ensure that the normal is of unit length and store it
  double len = sqrt(nx * nx + ny * ny + nz * nz);
  if(len > 1e-10) { 
    pair2Normal[currentSize][0] = nx / len;
    pair2Normal[currentSize][1] = ny / len;
    pair2Normal[currentSize][2] = nz / len;
  }
  else {
    filePrint(stderr, " ERROR: NodalContact::addCtcPair\n");
    filePrint(stderr, "   normal of ZERO length detected\n");
    exit(-1);
  }

  if(modePresent)
    if(Mode != 0) Mode = 1;
  if(fricCoefPresent) {
    pair2CoulombFrictionCoef[currentSize] = fricCoef;
    haveFriction = true;
    mode[currentSize] = 1;
  }
  else if(haveDefaultFriction) {
    haveFriction = true;
    if(modePresent) {
      mode[currentSize] = Mode;
      if(Mode == 0)
        pair2CoulombFrictionCoef[currentSize] = 1.0;
      else {
        if(globalMode == 1)
          pair2CoulombFrictionCoef[currentSize] = commonCoulombFrictionCoef;
        else
          pair2CoulombFrictionCoef[currentSize] = 0.0;
      }
    }
    else {
      mode[currentSize] = globalMode;
      pair2CoulombFrictionCoef[currentSize] = commonCoulombFrictionCoef;
    }
  }
  else {
    if(modePresent) {
      mode[currentSize] = Mode;
      if(Mode == 0)
        pair2CoulombFrictionCoef[currentSize] = 1.0;
      else
        pair2CoulombFrictionCoef[currentSize] = 0.0;
      // haveFriction = true;
    }
    else {
      pair2CoulombFrictionCoef[currentSize] = 0.0;
      mode[currentSize] = globalMode;
    }
  }
  gap0[currentSize]= normalGap;
  gap[currentSize] = normalGap;
  computeGap[currentSize] = !normalGapPresent;

  currentSize++;
 
  return 0;
}

//------------------------------------------------------
// Add virtual Elt in subToElem
int
NodalContact::addVirtualEltInSub(Connectivity **stn, int numSub) 
{
  // if (verboseFlag) std::cout << " ... Adding virtual elements to subdomains ..." << std::endl; 

  // Build nodeToSub
  Connectivity *&subToElem=*stn;
  Connectivity *elemToNode = new Connectivity(domain->getEset());
  Connectivity *subToNode = subToElem->transcon(elemToNode);
  Connectivity *nodeToSub = subToNode->reverse();
  
  // Build virtual nodes/elts and subToVirtualElt
  int *size = new int[numSub];  
  ResizeArray<int> **subToVirtualElt = new ResizeArray<int>*[numSub];
  int isub;
  for(isub = 0; isub < numSub; isub++) {
    size[isub] = 0;
    subToVirtualElt[isub] = new ResizeArray<int>(0);
  }

  int newNodeNb = (domain->getNodes()).last();
  if(iNodeStart < 0) iNodeStart = newNodeNb;
  double nothing[3] = {0}; // virtual coordinates;

  int newEltNb = (domain->getEset())->last();
  if(iEltStart < 0) iEltStart = newEltNb; // define this later, for the packed Elt set
  int eltNodesList[2];

  // GR:
  //  useMultipleLMP = false -> one and only one Lagrange multiplier is used
  //                            at each contact DOF
  //  useMultipleLMP = true  -> multiple Lagrange multipliers are used at
  //                            contact DOF's that describe contact between
  //                            points connected to more than one subdomain
  bool useMultipleLMP = true;
  int i;
  for(i = 0; i < pairsNb; i++) {
    int nL = pairsList[i][0], nR = pairsList[i][1];
    int numSubdL = 1;
    int numSubdR = 1;
    if(useMultipleLMP) {
      numSubdL = nodeToSub->num(nL);
      numSubdR = nodeToSub->num(nR);
    }    
    scaling[i] = 1.0e+00 / ( numSubdL + numSubdR);
    for(int subL = 0; subL < numSubdL; subL++) {
      int sdL = (*nodeToSub)[nL][subL];
      for(int subR = 0; subR < numSubdR; subR++) {
        int sdR = (*nodeToSub)[nR][subR];
        domain->addNode(newNodeNb,nothing);	  
        eltNodesList[0] = newNodeNb++; // the first node must be the virtual one
        eltNodesList[1] = nL; domain->addElem(newEltNb, virtualEltTypeNb, 2, eltNodesList); 
        (*(subToVirtualElt[sdL]))[size[sdL]++] = newEltNb++;
        eltNodesList[1] = nR; 
        domain->addElem(newEltNb, virtualEltTypeNb, 2, eltNodesList); 
        (*(subToVirtualElt[sdR]))[size[sdR]++] = newEltNb++;
      } // subR
    } // subL
  } // i
  
  if(newNodeNb != domain->getNodes().last() || newEltNb != domain->getEset()->last()) {
    std::cout << "Error in addVirtualEltInSub " << std::endl;
    exit(-1);
  }

  domain->setNumnodes(newNodeNb);

  // Build virtualNode2Pair  
  if(virtualNode2Pair) delete [] virtualNode2Pair;
  int numnodes = domain->getNodes().last();
  virtualNode2Pair = new int [numnodes];
  
  for(i=0; i<numnodes; i++) virtualNode2Pair[i] = -1;
  newNodeNb = iNodeStart;
  for(i=0; i<pairsNb; i++) {
    int nL = pairsList[i][0], nR = pairsList[i][1]; 
    int numSubdL = 1;
    int numSubdR = 1;
    if(useMultipleLMP) {
      numSubdL = nodeToSub->num(nL);
      numSubdR = nodeToSub->num(nR);
    }     
    for(int subL = 0; subL < numSubdL; subL++) {
      int sdL = (*nodeToSub)[nL][subL];
      for(int subR = 0; subR < numSubdR; subR++) {
        int sdR = (*nodeToSub)[nR][subR];
	virtualNode2Pair[newNodeNb++] = i;
      }
    }
  } // i
  
  // Build cx and connect for subToElem
  int *cx = new int[numSub+1];
  int *connect = new int[newEltNb];
  int curEle = 0,iele;
  for(isub=0; isub < numSub; ++isub) {
    cx[isub] = curEle;
    for(iele = 0; iele < subToElem->num(isub); ++iele) 
      connect[curEle++] = (*subToElem)[isub][iele];
    for(iele = 0; iele < size[isub]; ++iele)
      connect[curEle++] = (*(subToVirtualElt[isub]))[iele];    
  }
  cx[numSub] = curEle;
  delete subToElem; 
  subToElem = new Connectivity(numSub, cx, connect);

  // clean up
  if(elemToNode) delete elemToNode;
  if(subToNode) delete subToNode;
  if(nodeToSub) delete nodeToSub;
  if(subToVirtualElt) {
    for(i=0; i<numSub; ++i) if(subToVirtualElt[i]) delete subToVirtualElt[i];
    delete [] subToVirtualElt;
  }
  if(size) delete [] size;

  return 0;
}

//------------------------------------------------------
// Remove virtual Elt in subToElem
int
NodalContact::removeVirtualEltInSub(Connectivity **stn,int numSub) {

//  if (verboseFlag) 
//    std::cout << " ... Removing virtual contact elements from subdomains ..." << std::endl; 

  Connectivity *&subToElem=*stn;

  // Build cx and connect for subToElem
  int *cx = new int[numSub+1];
  int *connect = new int[iEltStart];
  int curEle = 0,iele;
  for(int isub=0; isub < numSub; ++isub) {
    cx[isub] = curEle;
    for(iele = 0; iele < subToElem->num(isub); ++iele) 
      if((*subToElem)[isub][iele] < iEltStart) 
	connect[curEle++] = (*subToElem)[isub][iele];
  }
  cx[numSub] = curEle;
  delete subToElem; 
  subToElem = new Connectivity(numSub, cx, connect);

  return 0;
}

//------------------------------------------------------
int
NodalContact::renum(int *table, int glmax, int numnodes, int *glNums) 
{
  // create the list of ctc couples in each sub
  // first : real node  (local sub nb) | second : sign
  int subPairsNb = 0;
  int *oldPairsList2New = new int[pairsNb];
  int (*SubPairsList)[2] = new int[pairsNb][2];

  for(int k=0; k<pairsNb; k++) {
    for(int j=0; j<2; j++) {
      if(pairsList[k][j] > glmax) SubPairsList[k][j] = -1;
      else SubPairsList[k][j] = table[pairsList[k][j]];
    }

    // the real node is now in pairsList[k][0] 
    // the sign involved in inequality [u(M)-u(M')].n <= gap is stored in pairsList[k][1]
    if (SubPairsList[k][0] >= 0) 
      SubPairsList[k][1] = 1;   // +u(M) ; = -1 if k not in local sub
    else { 
      SubPairsList[k][0] = SubPairsList[k][1];
      SubPairsList[k][1] = -1;  // -u(M')
    }

    if(SubPairsList[k][0] >= 0) oldPairsList2New[k] = subPairsNb++;
 
  } // k
  pairsList = SubPairsList; // don't delete : was a domain pointer

  bool *SubIsCtcNode = new bool[numnodes];
  int *SubVirtualNode2Pair = new int[numnodes];
  int i;
  for(i = 0; i < numnodes; i++) {
    if(glNums[i] < globalNumNodes)
      SubIsCtcNode[i] = isCtcNode[glNums[i]]; 
    else SubIsCtcNode[i] = false;
    if((glNums[i] >= domain->getLastNodeNb()) || (virtualNode2Pair[glNums[i]] < 0))  
      SubVirtualNode2Pair[i] = -1;
    else 
      SubVirtualNode2Pair[i] = oldPairsList2New[virtualNode2Pair[glNums[i]]];
  }
  isCtcNode = SubIsCtcNode; // don't delete : was a domain pointer
  virtualNode2Pair = SubVirtualNode2Pair;

  // count the number of duplicate contact pair numbers
  bool *process = new bool[numnodes];
  int numDuplicate = 0, l;
  for(i = 0; i < numnodes; i++) {
    if(virtualNode2Pair[i] < 0)
      process[i] = false;
    else
      process[i] = true;
  }
  for(i = 0; i < ( numnodes - 1); i++) {
    if(!process[i]) continue;
    for(l = i + 1; l < numnodes; l++) {
      if(virtualNode2Pair[l] == virtualNode2Pair[i]) {
        numDuplicate++;
        process[l] = false;
      }
    }
  }
  delete [] process;
  
  // set the number of contact pairs
  int numCtcDOFs = subPairsNb + numDuplicate;

  // keep only the contact nodes of the local sub
  resize<double,3> (pair2Normal, numCtcDOFs, subPairsNb);
  resize<double> (gap0, numCtcDOFs, subPairsNb);
  resize<double> (gap, numCtcDOFs, subPairsNb);
  resize<double> (scaling, numCtcDOFs, subPairsNb);
  resize<double> (pair2CoulombFrictionCoef, numCtcDOFs, subPairsNb);
  resize<bool> (computeGap, numCtcDOFs, subPairsNb);
  resize<int> (mode, numCtcDOFs, subPairsNb);
  pair2SlidingStatus = new int[numCtcDOFs][3];
  int kPair;
  for(kPair = 0; kPair < subPairsNb; kPair++) {
    this->setSlidingStatus(kPair, 1, 0);
  }
  if(haveFriction) {
    resize<double,3> (pair2Tangent1, numCtcDOFs, subPairsNb);
    resize<double,3> (pair2Tangent2, numCtcDOFs, subPairsNb);
    resize<double> (gapTg1, numCtcDOFs, subPairsNb);
    resize<double> (gapTg2, numCtcDOFs, subPairsNb);
    pair2TgExists = new bool[numCtcDOFs][2];
    pair2Tc1 = new double[numCtcDOFs][3];
    pair2Tc2 = new double[numCtcDOFs][3];
    pair2Tc1InActiveSet = new bool[numCtcDOFs];
    pair2Tc1NullifyPrimal = new bool[numCtcDOFs];
    pair2Tc2NullifyPrimal = new bool[numCtcDOFs];
    pair2Rconverged = new bool[numCtcDOFs];
    pair2SlipValue = new double[numCtcDOFs];
    int iPair;
    for(iPair = 0; iPair < subPairsNb; iPair++) {
      if(fabs(this->getPair2CoulombFricCoef(iPair)) < zeroForceTol) {
        this->setPair2Tc1InActiveSet(iPair, false);
        this->setPair2Tc1NullifyPrimal(iPair, true);
        this->setPair2Tc2NullifyPrimal(iPair, true);
        this->setSlidingStatus(iPair, 2, this->getDirChangeRange());
      }
      else {
        this->setPair2Tc1InActiveSet(iPair, true);
        this->setPair2Tc1NullifyPrimal(iPair, false);
        this->setPair2Tc2NullifyPrimal(iPair, false);
        this->setSlidingStatus(iPair, 2, (this->getDirChangeRange() + 1));
      }
      this->setSlipValue(iPair, this->getPair2CoulombFricCoef(iPair));
      this->setPair2Rconverged(iPair, false);
      this->setSlidingStatus(iPair, 0, -1);
    }
  }
  resize<int,2> (pairsList, numCtcDOFs, subPairsNb); // must be the last resize op
  pair2ccDofs = new int[numCtcDOFs][3];
  pair2cDofs = new int[numCtcDOFs][3];
  deltaFctc = new double [numCtcDOFs];
  // make room for the friction directions of deltaFctc
  if(haveFriction) {
    deltaFctcTg1 = new double [numCtcDOFs];
    deltaFctcTg2 = new double [numCtcDOFs];
  }
  pairsNb = numCtcDOFs;
  getSPatSize() = pairsNb;
  nullifyResidual = new bool[pairsNb];
  this->resetNullifyResidual();
  numConnectedNodes = new unsigned int[pairsNb];
  for(int iPair = 0; iPair < pairsNb; iPair++)
    numConnectedNodes[iPair] = 0;

  // GR:
  //  At this point, some contact Lagrange multipliers share the same 
  //  contact pairNb number. This is inconsistent with the setup
  //  used when running FetiDP-C from Salinas. Therefore, the
  //  next routine will ensure that each contact Lagrange multiplier
  //  will get its own (non-shared) pairNb number.
  createNonSharedCtcNumbers(subPairsNb, numnodes);

  // Construct local rotation matrices for scaling
  buildRotationMatrices();

  //std::cerr << "initializing globCtcIdentLists, pairsNb = " << pairsNb << std::endl;
  // globCtcIdentLists = new List<int> * [pairsNb];
  // for(i=0; i<pairsNb; ++i) globCtcIdentLists[i] = 0;

  localNumNodes = numnodes;  

// PJSA: DEBUG
/*
  if(pairsNb > 0) {
     cout << "SubDomain renumbered ctc " << endl;

     cout << "pairsNb = " << pairsNb << endl;
     cout << "pairsList (local numbering) = ";
     for(i = 0; i < pairsNb; i++)
       cout << pairsList[i][0] << " ";
     cout << endl;
     cout << "pairsList (global numbering) = ";
     for(i = 0; i < pairsNb; i++)
       cout << glNums[pairsList[i][0]] << " ";
     cout << endl;
     cout << "isCtcNode true (global numbering): ";
     int cc = 0;
     for(i = 0; i<localNumNodes; ++i) 
       if(isCtcNode[i]) { cout << glNums[i] << " "; cc++; }
     cout << endl;
     cerr << "numnodes = " << numnodes << ", glmax = " << glmax 
          << ", cc = " << cc << endl;
     cerr << "didn't find ";
     for(i = 0; i<nnodes; ++i)
       if(isCtcNode[i]) {
         bool found = false;
         for(int j = 0; j < pairsNb; j++)
           if(glNums[i] == glNums[pairsList[j][0]]) found = true;
         if(found == false) {
           cerr << glNums[i] << " ";
         }
       }
     cerr << " in pairsList \n";
  }
*/
// END DEBUG

  // clean up
  delete [] oldPairsList2New;
  delete [] SubPairsList;

  return 0;
}

//------------------------------------------------------
template <class Type, int d>
int
NodalContact::resize(Type (*&tab)[d], int dim, int checkDim) 
{
  // if(dim == 0) return 1;
  Type (*newTab)[d] = new Type[dim][d]; 
  if(dim == 0) { tab = newTab; return 1; }
  int size = 0;
  for(int i=0; i<pairsNb; i++) 
    if(pairsList[i][0] >= 0) {
      for(int j=0; j<d; j++) newTab[size][j] = tab[i][j];
      size++;
    }
  if(size != checkDim) { std::cout << "Error in NodalContact::resize : " << size << " " << checkDim << std::endl;exit(-1); }
  tab = newTab; // don't delete : was a domain pointer

  return 0;
}

//------------------------------------------------------
template <class Type>
int
NodalContact::resize(Type *&tab, int dim, int checkDim) 
{
  // if(dim == 0) return 1;
  Type *newTab = new Type[dim];   
  if(dim == 0) { tab = newTab; return 1; }
  int size = 0;
  for(int i=0; i<pairsNb; i++) 
    if(pairsList[i][0] >= 0) 
      newTab[size++] = tab[i];
  if(size != checkDim) { std::cout << "Error in NodalContact::resize : " << size << " " << checkDim << std::endl;exit(-1); }
  tab = newTab; // don't delete : was a domain pointer
  return 0;
}

//------------------------------------------------------
int 
NodalContact::computeInitialGap(CoordSet &nodes) 
{
  if(haveFriction) {
    gapTg1 = new double[pairsNb];
    gapTg2 = new double[pairsNb];
  }
  int i;
  for (i=0;i<pairsNb;i++) {
    Node &n1=*nodes[pairsList[i][0]];
    Node &n2=*nodes[pairsList[i][1]];
    if(computeGap[i]) {
      gap0[i]= (n2.x-n1.x)*pair2Normal[i][0]+
               (n2.y-n1.y)*pair2Normal[i][1]+
               (n2.z-n1.z)*pair2Normal[i][2];
      gap[i] = gap0[i];
    }
    if(haveFriction) {
      gapTg1[i] = 0.0;
      gapTg2[i] = 0.0;
    }
  } // i   

  // NOTE: if one or more dofs involved in the contact have a non-zero dirichlet 
  // boundary condition then the gap should be modified accordingly
  // see DecDomain::distributeMPCs() for similar treatment of MPC rhs

  return 0;
}

//------------------------------------------------------
int
NodalContact::computeCtcPrimalError(int *allBoundDofs, int interfSize, double *u, int *boundDofFlag) 
{
  ctcPrimalError = 0;
  int iDof;
  for(iDof = 0; iDof < interfSize; ++iDof)
    if(boundDofFlag[iDof] == 1)   // PJSA: modified to ignore MPCs
      if(fabs(u[iDof]) > ctcPrimalError) ctcPrimalError = fabs(u[iDof]);
  
  return 0;
}

//------------------------------------------------------
void 
NodalContact::setInvBoundMap(int len, int *tab) 
{
  if(invBoundMap) delete[] invBoundMap;
  invBoundMap = new int[len];
  int i;
  for(i=0; i<len; i++) invBoundMap[i] = tab[i];
}

//------------------------------------------------------
void 
NodalContact::setInvInternalMap(int len,int * tab) 
{
  if(invInternalMap) delete[] invInternalMap;
  invInternalMap = new int[len];
  for(int i=0; i<len; i++) invInternalMap[i] = tab[i];
}

//------------------------------------------------------
void
NodalContact::initInactiveSetResidual() 
{
  if(inactiveSetResidual == 0)
    inactiveSetResidual = new ResizeArray<double>(0);
  inactiveSetSize=0;
}

//------------------------------------------------------
double
NodalContact::saveAndCheckInactiveSetResidual(double value) 
{
  double return_value = (*inactiveSetResidual)[inactiveSetSize];
  (*inactiveSetResidual)[inactiveSetSize++] = value;
  return return_value;
}

//------------------------------------------------------
double
NodalContact::getInactiveSetResidual() 
{
  if(inactiveSetSize)
    return (*inactiveSetResidual)[--inactiveSetSize]; 
  else {
    std::cout << "Error in getInactiveSetResidual " << std::endl;
    exit(-1); return 0.0;
  }
}

//------------------------------------------------------
void
NodalContact::initInactiveSetFrictionResidual()
{
  if(inactiveSetFrictionResidual == 0) 
    inactiveSetFrictionResidual = new ResizeArray<double>(0);
  inactiveSetFrictionSize = 0;
}

//------------------------------------------------------
void
NodalContact::saveInactiveSetFrictionResidual(double Tg1Value, 
                                              double Tg2Value)
{
  ( *inactiveSetFrictionResidual)[inactiveSetFrictionSize++] = Tg1Value;
  ( *inactiveSetFrictionResidual)[inactiveSetFrictionSize++] = Tg2Value;
}

//------------------------------------------------------
double
NodalContact::getInactiveSetFrictionResidual()
{
  if(inactiveSetFrictionSize)
    return (*inactiveSetFrictionResidual)[--inactiveSetFrictionSize]; 
  else {
    std::cout << "Error in getInactiveSetFrictionResidual " << std::endl;
    exit(-1);
    return 0.0;
  }
}

//------------------------------------------------------
void
NodalContact::setDof2Pair(int nodeNb, int dofNb, int nnn) {

  if(dof2Pair==0) dof2Pair= new ResizeArray<int>(-1);
  int j = -1;
  int i;
  for(i = 0; i < pairsNb; i++)
    if(pairsList[i][0] == nodeNb) {j=i; break;}
  if(j < 0) {
     std::cout << "Error in setDof2Pair " << std::endl;
     std::cout << "Have " << pairsNb << " pairs, looking for " << nodeNb << std::endl;
     for(int i = 0; i < pairsNb; i++)
       std::cout << pairsList[i][0] << " ";
     std::cout << std::endl;
     exit(-1);
  }
  for (int k = 0; k < nnn; k++) (*dof2Pair)[dofNb+k] = j;
}

//------------------------------------------------------
int 
NodalContact::getDof2Pair(int dofNb) 
{
  if(dof2Pair == 0) return -1;
  return (*dof2Pair)[dofNb];
}

//------------------------------------------------------
int 
NodalContact::buildRotationMatrices() 
{
  // necessary for scaling (see multKbbCtc in Subdomain class)
  // n is the normal
  // Compute the tangents t1 and t2 = n ^ t1
  // The rotation matrice R is defined by transpose(R) = |transpose(n)          | 
  //                                                     |coeff * transpose(t1) | 
  //                                                     |coeff * transpose(t2) |
  // R e1 = n
  // R e2 = t1
  // R e3 = t2
  // R = | R0 R1 R2 |
  //     | R3 R4 R5 |
  //     | R6 R7 R8 |
  const double epsilon = 1E-6;
  double *n,*rot;
  double t1[3],t2[3],coeff;

  if(pair2rotMat) delete [] pair2rotMat;
  pair2rotMat = new double [pairsNb][9];

  for(int i = 0; i < pairsNb; i++) {
    n = pair2Normal[i];
    if((1 - n[2] * n[2]) > epsilon) { //  -> can divide by 1-n[2]*n[2]
      t1[0] = -n[1]; t1[1] = n[0]; t1[2] = 0;
      t2[0] = -n[0] * n[2]; t2[1] = -n[1] * n[2]; t2[2] = 1 - n[2] * n[2];
      coeff = 1 / (1 - n[2] * n[2]);      
    }
    else {                       // -> can divide by 1-n[1]*n[1]
      t1[0]=n[2];t1[1]=0;t1[2]=-n[0]; // not tested
      t2[0]=-n[0]*n[1];t2[1]=1-n[1]*n[1];t2[2]=-n[1]*n[2];
      coeff=1/(1-n[1]*n[1]);
    }
    
    rot=pair2rotMat[i];
    rot[0]=n[0];rot[3]=n[1];rot[6]=n[2];
    rot[1]=coeff*t1[0];rot[4]=coeff*t1[1];rot[7]=coeff*t1[2];
    rot[2]=coeff*t2[0];rot[5]=coeff*t2[1];rot[8]=coeff*t2[2];

  } // i

  return 0;

}

//------------------------------------------------------
int
NodalContact::rotateTranspose(double *u,double *tRu,int i) {
  double *rot = pair2rotMat[i];
  tRu[0]=rot[0]*u[0]+rot[3]*u[1]+rot[6]*u[2];
  tRu[1]=rot[1]*u[0]+rot[4]*u[1]+rot[7]*u[2];
  tRu[2]=rot[2]*u[0]+rot[5]*u[1]+rot[8]*u[2];
  return 0;
}

//------------------------------------------------------
int
NodalContact::rotate(double *u,double *Ru,int i) {
  double *rot = pair2rotMat[i];
  Ru[0]=rot[0]*u[0]+rot[1]*u[1]+rot[2]*u[2];
  Ru[1]=rot[3]*u[0]+rot[4]*u[1]+rot[5]*u[2];
  Ru[2]=rot[6]*u[0]+rot[7]*u[1]+rot[8]*u[2];
  return 0;
}

//------------------------------------------------------
// Compute M = rotCenter + RotMat (M-rotCenter)
// rotation matrix in 0:8 ; center of rotation in 9:11
int
NodalContact::rotate(double **M) {
  double &mx=(*M)[0]; double &my=(*M)[1]; double &mz=(*M)[2];
  double &cx=rotMat[9];double &cy=rotMat[10];double &cz=rotMat[11];
  double *m=rotMat;
  double x=cx+m[0]*(mx-cx)+m[1]*(my-cy)+m[2]*(mz-cz);
  double y=cy+m[3]*(mx-cx)+m[4]*(my-cy)+m[5]*(mz-cz);
  double z=cz+m[6]*(mx-cx)+m[7]*(my-cy)+m[8]*(mz-cz);
  mx=x;my=y;mz=z;

  return 0;
}

//------------------------------------------------------
// Compute rotMat Vec
int
NodalContact::rotate(double *Mx,double *My,double *Mz) {
  double &mx=*Mx; double &my=*My; double &mz=*Mz;
  double *m=rotMat;
  double x=m[0]*mx+m[1]*my+m[2]*mz;
  double y=m[3]*mx+m[4]*my+m[5]*mz;
  double z=m[6]*mx+m[7]*my+m[8]*mz;
  mx=x;my=y;mz=z;
  
  return 0;
}

//------------------------------------------------------
void
NodalContact::nodeMerging(int num, int etype, int nnodes, int *n)
{
  int inode, jnode;
  for(inode = 0; inode < nnodes; inode++) {
    for(jnode = 0; jnode < pairsNb; jnode++)
      if(n[inode] == pairsList[jnode][1]) {
	// std::cout << "__CTC__ Merging " << pairsList[jnode][1] << " with " << pairsList[jnode][0] << std::endl;
	n[inode] = pairsList[jnode][0]; // always keep the first number	
	break;
      }
  }
}

//------------------------------------------------------
void
NodalContact::setUnitTangents()
{
  if(haveFriction) {
    int istat;
    unsigned int nLargestIndex, nSmallestIndex;
    double *normal, *tangent1, *tangent2, nLargestValue,
           nSmallestValue, e1[3], e2[3], e3[3], Q[9];
    pair2Tangent1 = new double[pairsNb][3];
    pair2Tangent2 = new double[pairsNb][3];
    for(int iPair = 0; iPair < pairsNb; iPair++) {
      normal = this->getPair2Normal( iPair);
      tangent1 = this->getPair2Tangent1( iPair);
      tangent2 = this->getPair2Tangent2( iPair);
      // if 2D contact -> choose e3 to normal transformation such that
      //                  e3 lies in the 2D problem plane
      // further ensure that dot(e3,normal) != -1.0
      nSmallestIndex = 0;
      nSmallestValue = fabs( normal[0]);
      int j;
      for(j = 1; j < 3; j++) {
        if(fabs( normal[j]) < nSmallestValue) {
          nSmallestValue = fabs(normal[j]);
          nSmallestIndex = j;
        }
      }
      nLargestIndex = 0;
      nLargestValue = -2.0;
      for(j = 0; j < 3; j++) {
        if(j != nSmallestIndex) {
          if(normal[j] > nLargestValue) {
            nLargestValue = normal[j];
            nLargestIndex = j;
          }
        }
      }
      switch(nLargestIndex) {
        case 0:
          e1[0] = 0.0; e1[1] = 1.0; e1[2] = 0.0;
          e2[0] = 0.0; e2[1] = 0.0; e2[2] = 1.0;
          e3[0] = 1.0; e3[1] = 0.0; e3[2] = 0.0;
          break;
        case 1:
          e1[0] = 0.0; e1[1] = 0.0; e1[2] = 1.0;
          e2[0] = 1.0; e2[1] = 0.0; e2[2] = 0.0;
          e3[0] = 0.0; e3[1] = 1.0; e3[2] = 0.0;
          break;
        case 2:
          e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0;
          e2[0] = 0.0; e2[1] = 1.0; e2[2] = 0.0;
          e3[0] = 0.0; e3[1] = 0.0; e3[2] = 1.0;
          break;
      }
      istat = b2QEtoA( e3, normal, Q);
      if(istat == -2) {
        filePrint(stderr, " ERROR: NodalContact::setUnitTangents\n");
        filePrint(stderr, "   division by zero detected while constructing");
        filePrint(stderr, " the zero drill rotation matrix\n");
        exit(-1);
      }
      b2MatrixProduct(3, 3, 1, Q, e1, tangent1);
      b2MatrixProduct(3, 3, 1, Q, e2, tangent2);
    }
  }
}

void
NodalContact::reference2FrictionBasis(int i, double ug, double vg, 
                                      double &uc, double &vc)
{
  double *tg1 = this->getPair2Tangent1(i);
  double *tc1 = this->getPair2Tc1(i);
  double *tc2 = this->getPair2Tc2(i);
  double normTg1 = sqrt( tg1[0] * tg1[0] + tg1[1] * tg1[1] + tg1[2] * tg1[2]);
  double normTc1 = sqrt( tc1[0] * tc1[0] + tc1[1] * tc1[1] + tc1[2] * tc1[2]);
  double arg = (tg1[0] * tc1[0] + tg1[1] * tc1[1] + tg1[2] * tc1[2]) /
               (normTg1 * normTc1);
  if(arg >  1.0) arg =  1.0;
  if(arg < -1.0) arg = -1.0;
  double alpha = acos(arg);
  if((tg1[0] * tc2[0] + tg1[1] * tc2[1] + tg1[2] * tc2[2]) > 0.0)
    alpha = -alpha;
  uc = cos(alpha) * ug + sin(alpha) * vg;
  vc = -sin(alpha) * ug + cos(alpha) * vg;
}

void
NodalContact::friction2ReferenceBasis(int i, double uc, double vc, 
                                      double &ug, double &vg)
{
  double *tg1 = this->getPair2Tangent1(i);
  double *tc1 = this->getPair2Tc1(i);
  double *tc2 = this->getPair2Tc2(i);
  double normTg1 = sqrt(tg1[0] * tg1[0] + tg1[1] * tg1[1] + tg1[2] * tg1[2]);
  double normTc1 = sqrt(tc1[0] * tc1[0] + tc1[1] * tc1[1] + tc1[2] * tc1[2]);
  double arg = (tg1[0] * tc1[0] + tg1[1] * tc1[1] + tg1[2] * tc1[2]) /
               (normTg1 * normTc1);
  if(arg >  1.0) arg =  1.0;
  if(arg < -1.0) arg = -1.0;
  double alpha = acos(arg);
  if((tg1[0] * tc2[0] + tg1[1] * tc2[1] + tg1[2] * tc2[2]) > 0.0)
    alpha = -alpha;
  ug = cos(alpha) * uc - sin(alpha) * vc;
  vg = sin(alpha) * uc + cos(alpha) * vc;
}

void 
NodalContact::constructTc2(int i)
{
  double *normal = this->getPair2Normal(i);
  double normalNorm = b2VectorNorm(3, &normal[0]);
  double rotationVector[3], work[9], q[9];
  double *tg1 = this->getPair2Tangent1(i);
  double *tg2 = this->getPair2Tangent2(i);
  double *tc1 = this->getPair2Tc1(i);
  double *tc2 = this->getPair2Tc2(i);
  double normTg1 = sqrt(tg1[0] * tg1[0] + tg1[1] * tg1[1] + tg1[2] * tg1[2]);
  double normTc1 = sqrt(tc1[0] * tc1[0] + tc1[1] * tc1[1] + tc1[2] * tc1[2]);
  double arg = (tg1[0] * tc1[0] + tg1[1] * tc1[1] + tg1[2] * tc1[2]) /
               (normTg1 * normTc1);
  if(arg >  1.0) arg =  1.0;
  if(arg < -1.0) arg = -1.0;
  double alpha = acos(arg);
  if((tc1[0] * tg2[0] + tc1[1] * tg2[1] + tc1[2] * tg2[2]) < 0.0)
    alpha = -alpha;
  int k;
  for(k = 0; k < 3; k++)
    rotationVector[k] = alpha * normal[k] / normalNorm;
  b2RotationMatrix(&rotationVector[0], &work[0], &q[0]);
  b2MatrixProduct(3, 3, 1, &q[0], &tg2[0], &tc2[0]);
}

void 
NodalContact::createNonSharedCtcNumbers(int sharedDim, int numNodes)
{
  // local objects
  int i, j, k, iPair, strNext = sharedDim;

  for(i = 0; i < numNodes - 1; i++) {
    iPair = virtualNode2Pair[i];
    if(iPair < 0) continue;
    for(j = i + 1; j < numNodes; j++) {
      if(virtualNode2Pair[j] == iPair) {
        if(strNext < pairsNb) {
          virtualNode2Pair[j] = strNext;
          pairsList[strNext][0] = pairsList[iPair][0];
          pairsList[strNext][1] = pairsList[iPair][1];
          pair2Normal[strNext][0] = pair2Normal[iPair][0];
          pair2Normal[strNext][1] = pair2Normal[iPair][1];
          pair2Normal[strNext][2] = pair2Normal[iPair][2];
          gap0[strNext]= gap0[iPair];
          gap[strNext] = gap[iPair];
          computeGap[strNext] = computeGap[iPair];
          mode[strNext] = mode[iPair];
          scaling[strNext] = scaling[iPair];
          pair2CoulombFrictionCoef[strNext] = pair2CoulombFrictionCoef[iPair];
          pair2SlidingStatus[strNext][1] = pair2SlidingStatus[iPair][1];
          if(hasFriction()) {
            for(k = 0; k < 3; k++) {
              pair2Tangent1[strNext][k] = pair2Tangent1[iPair][k];
              pair2Tangent2[strNext][k] = pair2Tangent2[iPair][k];
            }
            gapTg1[strNext] = gapTg1[iPair];
            gapTg2[strNext] = gapTg2[iPair];
            pair2Tc1InActiveSet[strNext] = pair2Tc1InActiveSet[iPair];
            pair2Tc1NullifyPrimal[strNext] = pair2Tc1NullifyPrimal[iPair];
            pair2Tc2NullifyPrimal[strNext] = pair2Tc2NullifyPrimal[iPair];
            pair2SlipValue[strNext] = pair2SlipValue[iPair];
            pair2Rconverged[strNext] = pair2Rconverged[iPair];
            pair2SlidingStatus[strNext][0] = pair2SlidingStatus[iPair][0];
            pair2SlidingStatus[strNext][2] = pair2SlidingStatus[iPair][2];
          }
          strNext++;
        }
        else {
          std::cout << "Error in NodalContact::createNonSharedCtcNumbers : "
               << " more local contact pairs then expected"
               << std::endl;
          exit(-1);
        }
      }
    }
  }
}

void
NodalContact::initPlanCommSetupPointer()
{
  planCommSetup = 0;
}

PlanCommSetup &
NodalContact::getPlanCommSetup(int pairNb)
{
  return planCommSetup[pairNb];
}

void
NodalContact::allocatePlanCommSetup()
{
  if(planCommSetup != 0) delete [] planCommSetup;
  planCommSetup = new PlanCommSetup[pairsNb];
}

void 
NodalContact::setPlanCommSetupGlobStrNodes(int *glNums, int numnodes)
{
  int i, pairNb;
  if(pairsNb > 0) {
    for(i = 0; i < numnodes; i++) {
      pairNb = getVirtualNodeNb(i);
      if(pairNb >= pairsNb) std::cerr << "!!!error: pairNb = " << pairNb << ", pairsNb = " << pairsNb << std::endl;
      if(pairNb >= 0) 
        planCommSetup[pairNb].setGlobStrNode(glNums[pairsList[pairNb][0]]);
    }
  }
}

int &
NodalContact::getSPatSize()
{
  return sPatSize;
}
/*
// GR: original version
void
NodalContact::findGlobCtcIdent(int globStrNode, int globEndNode,
                               list<int> &globCtcIdentList, int numnodes)
{
  int i, l, pairNb, pairNb2, cnt1, cnt2;
  cnt1 = 0;
  for(i = 0; i < numnodes; i++) {
    pairNb = getVirtualNodeNb(i);
    cnt1++;
    if(pairNb >= 0) {
      cnt2 = 0;
      for(l = 0; l < numnodes; l++) {
        pairNb2 = getVirtualNodeNb(l);
        cnt2++;
        if((cnt2 < cnt1) && (pairNb == pairNb2)) break;
        if(cnt2 >= cnt1) break;
      }
      if((cnt2 < cnt1) && (pairNb == pairNb2)) continue;
      if((globStrNode == planCommSetup[pairNb].getGlobStrNode()) &&
         (globEndNode == planCommSetup[pairNb].getGlobEndNode())) {
        globCtcIdentList.push_back(pairNb + 1); 
        continue;
      }
      if((globEndNode == planCommSetup[pairNb].getGlobStrNode()) &&
         (globStrNode == planCommSetup[pairNb].getGlobEndNode())) {
        globCtcIdentList.push_back(-pairNb - 1);
      }
    }
  }
}
*/
// PJSA: new version
void
NodalContact::findGlobCtcIdent(int globStrNode, int globEndNode,
                               list<int> &globCtcIdentList, int printFlag)
{
  int i;
  globCtcIdentList.clear();
  for(i = 0; i < pairsNb; i++) {
    if((globStrNode == planCommSetup[i].getGlobStrNode()) &&
       (globEndNode == planCommSetup[i].getGlobEndNode())) {
      globCtcIdentList.push_back(i + 1);  
    }
    else if((globEndNode == planCommSetup[i].getGlobStrNode()) &&
            (globStrNode == planCommSetup[i].getGlobEndNode())) {
      globCtcIdentList.push_back(-i - 1);
    }
  }

 if(printFlag == -1) {
  std::cerr << "globStrNode = " << globStrNode << ", globEndNode = " << globEndNode
       << ", globCtcIdentList = ";
  list<int>::iterator listItr;
  for(listItr = globCtcIdentList.begin(); listItr != globCtcIdentList.end(); ++listItr)
    std::cerr << *listItr << " ";
  std::cerr << std::endl;
 }
}

void NodalContact::getGlobCtcIdent(int globStrNode, int globEndNode, list<int> &fndList, int numnodes)
{
  CtcMapType::iterator iter = globCtcIdentMap.find(pair<int, int>(globStrNode, globEndNode));
  if(iter != globCtcIdentMap.end()) {
    // std::cerr << "found list \n";
    // std::cerr << "globStrNode = " << globStrNode << ", globEndNode = " << globEndNode << std::endl;
    fndList = (*iter).second;
  }
  else {
    // std::cerr << "making list \n";
    findGlobCtcIdent(globStrNode, globEndNode, fndList, numnodes);
    globCtcIdentMap.insert(CtcValuePair(pair<int, int>(globStrNode, globEndNode), fndList));
  }
}

void
PlanCommSetup::showData()
{
  std::cerr << "globStrNode = " << globStrNode << ", globEndNode = " << globEndNode 
       << ", communicate = " << communicate << ", subdList = ";
  list<int>::iterator listItr;
  for(listItr = subdList.begin(); listItr != subdList.end(); ++listItr)  
    std::cerr << *listItr << " ";
  std::cerr << std::endl;
}

double*
NodalContact::getNormalOrTangent(int pairNb, int pairTg)
{
  double *ret;
  switch(pairTg) {
    default:
    case 0:
      ret = getPair2Normal(pairNb);
      break;
    case 1:
      ret = (getPair2TgExists(pairNb)[0]) ? getPair2Tangent1(pairNb) : getPair2Tangent2(pairNb);
      break;
    case 2:
      ret = getPair2Tangent2(pairNb);
      break;
  }
  return ret;
}

void 
NodalContact::print(int* subNum)
{
  if(subNum)
       fprintf(stderr," - NodalContact object of subd %3d, contains %4d contact pairs:\n",*subNum, pairsNb); 
  else fprintf(stderr," - NodalContact object contains %4d contact pairs:\n",pairsNb); 
  fflush(stderr);
  for(int i=0; i<pairsNb; i++) {
    fprintf(stderr,"   * pair %4d, (local) node %6d, normal = %2e %2e %2e, gap = %2e, sign = %d\n",
            i, pairsList[i][0],pair2Normal[i][0],pair2Normal[i][1],pair2Normal[i][2],gap[i],int(dsign(i)));
    fflush(stderr);
  }
}

int * 
NodalContact::operator[] (int i) { return pairsList[i]; };

bool 
NodalContact::getFlagCtc() { return ((pairsList != 0) || flagCtc); }; 

void
NodalContact::setFlagCtc(bool flag) { flagCtc = flag; };

int 
NodalContact::getPairsNb() { return pairsNb; };

int
NodalContact::getPairsNbTotal() { return pairsNbTotal; };

int 
NodalContact::virtualElt(int i,int j) { return iEltStart+2*i+j; };  

double * 
NodalContact::getPair2Normal(int i) { return pair2Normal[i]; };


double * 
NodalContact::getPair2Tangent1(int i)
{
  if(haveFriction)
    return pair2Tangent1[i];
  else
    return 0;
}

double * 
NodalContact::getPair2Tangent2(int i)
{
  if(haveFriction)
    return pair2Tangent2[i];
  else
    return 0;
}

double * 
NodalContact::getPair2Tc1(int i)
{
  if(haveFriction)
    return pair2Tc1[i];
  else
    return 0;

}

double * 
NodalContact::getPair2Tc2(int i)
{
  if(haveFriction)
    return pair2Tc2[i];
  else
    return 0;
}

bool *
NodalContact::getPair2TgExists(int i)
{
  if(haveFriction)
    return pair2TgExists[i];
  else
    return 0;
}

void 
NodalContact::setPair2Tc1InActiveSet(int i, bool active)
{
  if(haveFriction) pair2Tc1InActiveSet[i] = active;
}
  
void
NodalContact::setPair2Tc1NullifyPrimal(int i, bool nullify)
{
  if(haveFriction) pair2Tc1NullifyPrimal[i] = nullify;
}

void
NodalContact::setPair2Tc2NullifyPrimal(int i, bool nullify)
{
  if(haveFriction) pair2Tc2NullifyPrimal[i] = nullify;
}

void 
NodalContact::setPair2Rconverged(int i, bool converged)
{
  if(haveFriction) pair2Rconverged[i] = converged;
}

bool
NodalContact::isTc1InActiveSet(int i)
{
  if(haveFriction)
    return pair2Tc1InActiveSet[i];
  else
    return false;
}

bool
NodalContact::isRconverged(int i)
{
  if(haveFriction)
    return pair2Rconverged[i];
  else
    return false;
}

bool 
NodalContact::isTc1NullifiedInPrimal(int i)
{
  if(haveFriction)
    return pair2Tc1NullifyPrimal[i];
  else
    return false;
}

bool 
NodalContact::isTc2NullifiedInPrimal(int i)
{
  if(haveFriction)
    return pair2Tc2NullifyPrimal[i];
  else
    return false;
}

bool 
NodalContact::IsVirtualCtcNode(int NodeNb) {return ((virtualNode2Pair!=0) && (pairsNb!=0) &&
						    (virtualNode2Pair[NodeNb]>=0));} 
bool 
NodalContact::IsCtcNode(int inode) {return  ((isCtcNode!=0) && (pairsNb!=0) && isCtcNode[inode]); };

int 
NodalContact::getVirtualNodeNb(int NodeNb) {return virtualNode2Pair[NodeNb];};

int *
NodalContact::getPair2ccDofs(int pairNb) {return pair2ccDofs[pairNb];};

int *
NodalContact::getPair2cDofs(int pairNb) {return pair2cDofs[pairNb];};

int
NodalContact::setPair2Dofs(int i,int *dofccNb,int * dofcNb) {
  pair2ccDofs[i][0]=dofccNb[0];pair2ccDofs[i][1]=dofccNb[1];pair2ccDofs[i][2]=dofccNb[2]; 
  pair2cDofs[i][0]=dofcNb[0];pair2cDofs[i][1]=dofcNb[1];pair2cDofs[i][2]=dofcNb[2]; 
  return 0;}; 

void 
NodalContact::changeSign(int pair) {pairsList[pair][1] = -pairsList[pair][1];}

int
NodalContact::sign(int pair) {if (pairsList[pair][1]>=0) return 1; else return -1;};

double
NodalContact::dsign(int pair) {if (pairsList[pair][1]>=0) return 1.0; else return -1.0;};

double &
NodalContact::Gap(int pairNb, int pairTg)
{
  double *Gap;
  switch( pairTg) {
    case 0:
      Gap = gap;
      break;
    case 1:
      Gap = gapTg1;
      break;
    case 2:
      Gap = gapTg2;
      break;
  }
  return Gap[pairNb];
}

double
NodalContact::setNormalGap(int pairNb, double ngap)
{
  gap0[pairNb] = ngap;
  gap[pairNb]  = ngap;
  return(gap[pairNb]);
}

double 
NodalContact::updateNormalGap(int pairNb, double Dg, bool first)
{
  if(first) gap[pairNb] = gap0[pairNb] - Dg;
  else gap[pairNb] -= Dg;
  return(gap[pairNb]);
}

void 
NodalContact::setNormal(int pairNb, double* nrm)
{
  //ensure that the normal is of unit length and store it
  double len = sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2]);
  if(len > 1e-10) {
    pair2Normal[pairNb][0] = nrm[0]/ len;
    pair2Normal[pairNb][1] = nrm[1]/ len;
    pair2Normal[pairNb][2] = nrm[2]/ len;
  }
  else {
    cerr << " *** WARNING: normal of ZERO length detected\n";
  }
}

bool 
NodalContact::IsInActiveSet(int pair) 
{
  if(pair < 0) return false; 
  return (pairsList[pair][0] >= 0.0);
}

int 
NodalContact::changeStatus(int pair) 
{
  pairsList[pair][0] = -1 - pairsList[pair][0];
  return 0;
}

double
NodalContact::getCtcPrimalError() { return ctcPrimalError; }

double 
NodalContact::getEpsilonPrimal() { return epsilonPrimal; }

double
NodalContact::getEpsilonDual() { return epsilonDual; }

double
NodalContact::getCtcError0() 
{
  if(ctcError0 < 0) {
    std::cout << "CtcError0 max not initialized" << endl;
    exit(-1);
  } 
  return ctcError0;
}

void
NodalContact::setCtcError0(double ctcError) { ctcError0 = ctcError; }

int*
NodalContact::getInvBoundMap() { return invBoundMap; }

int*
NodalContact::getInvInternalMap() { return invInternalMap; }

void
NodalContact::setCtcPrimalError(double error) { ctcPrimalError = error; }

double
NodalContact::getPair2Scaling(int pairNb) { return scaling[pairNb]; }

bool
NodalContact::IsCtcYes(int pairNb) { return true; }

double
NodalContact::getFetiCtcTol() { return fetiCtcTol; }

void
NodalContact::setDeltaFctc(int pairNb, int pairTg, double val) 
{
  switch(pairTg) {
    case 0:
      deltaFctc[pairNb] = val;
      break;
    case 1:
      deltaFctcTg1[pairNb] = val;
      break;
    case 2:
      deltaFctcTg2[pairNb] = val;
      break;
  }
}

double
NodalContact::getDeltaFctc(int pairNb, int pairTg) 
{
  double val = 0.0;
  switch(pairTg) {
    case 0:
      val = deltaFctc[pairNb];
      break;
    case 1:
      val = deltaFctcTg1[pairNb];
      break;
    case 2:
      val = deltaFctcTg2[pairNb];
      break;
  }
  return val;
}

// Preprocessing
void 
NodalContact::setNodeMerging(int pair) { mergingNodes = true; setPairsNb(pair); }

bool 
NodalContact::nodeMerging() { return ((pairsList != 0) && mergingNodes); } 

int
NodalContact::setRotMat(double d0,double d1,double d2,
			double d3,double d4,double d5,
			double d6,double d7,double d8, // rotation matrix
			double d9,double d10,double d11) // center of rotation
{
  if (rotMat) delete [] rotMat; 
  rotMat = new double[12];
  rotMat[0] = d0; rotMat[1] = d1; rotMat[2] = d2; rotMat[3] = d3; 
  rotMat[4] = d4; rotMat[5] = d5; rotMat[6] = d6; rotMat[7] = d7;
  rotMat[8] = d8; rotMat[9] = d9; rotMat[10] = d10; rotMat[11] = d11;
  return 0;
}

bool
NodalContact::rotationDefined() { return rotMat != 0; }

int
NodalContact::setPairsNbTotal(int n) { pairsNbTotal = n; return 0; }

void 
NodalContact::setNullifyResidual(int pair, bool setValue) { nullifyResidual[pair] = setValue; }

bool 
NodalContact::isNullifiedInResidual(int pair) { return nullifyResidual[pair]; }

bool
NodalContact::hasFriction() { return haveFriction; }

double &
NodalContact::getPair2CoulombFricCoef(int pairNb) { return pair2CoulombFrictionCoef[pairNb]; }

double
NodalContact::getZeroForceTol() {return zeroForceTol; }

unsigned int &
NodalContact::getNumConnectedNodes(int pairNb) { return numConnectedNodes[pairNb]; }

void
NodalContact::setSlidingStatus(int i, int component, int setValue)
{ pair2SlidingStatus[i][component] = setValue; }

int 
NodalContact::getSlidingStatus(int i, int component)
{ return pair2SlidingStatus[i][component]; }

void 
NodalContact::setSlipValue(int i, double slipValue)
{ pair2SlipValue[i] = slipValue; }

double
NodalContact::getSlipValue(int i) { return pair2SlipValue[i]; }

int
NodalContact::getDirChangeRange()
{
  // the value returned by this function must be an odd and positive number,
  // i.e, 1,3,5,7...
  // It sets a range that is used to change the direction of the frictional
  // force in steps in case this force is not opposite to the slip
  return 5;
}

int
NodalContact::getMode(int pairNb) { return mode[pairNb]; }

void
PlanCommSetup::setGlobStrNode(int _globStrNode) { globStrNode = _globStrNode; }

void
PlanCommSetup::setGlobEndNode(int _globEndNode) {globEndNode = _globEndNode; }

int
PlanCommSetup::getGlobStrNode() { return globStrNode; }

int
PlanCommSetup::getGlobEndNode() { return globEndNode; }

list<int> &
PlanCommSetup::getSubdList() { return subdList; }

void
PlanCommSetup::setCommunicate(bool _communicate) { communicate = _communicate; }

bool
PlanCommSetup::getCommunicate() { return communicate; }

void
PlanCommSetup::initValue()
{
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
}

void
PlanCommSetup::addToValue(int component, double addValue)
{ value[component] += addValue; }

double *
PlanCommSetup::getValue() { return value; }

void
PlanCommSetup::initTotalValue()
{
  totalValue[0] = 0.0;
  totalValue[1] = 0.0;
  totalValue[2] = 0.0;
}

void
PlanCommSetup::addToTotalValue(double *addValue)
{
  totalValue[0] += addValue[0];
  totalValue[1] += addValue[1];
  totalValue[2] += addValue[2];
}

double
PlanCommSetup::getTotalValue(int component) { return totalValue[component]; }

void
PlanCommSetup::initMultiplicity() { multiplicity = 0; }

void
PlanCommSetup::addToMultiplicity(int addValue) { multiplicity += addValue; }

int
PlanCommSetup::getMultiplicity() { return multiplicity; }

