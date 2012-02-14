#include <Element.d/SuperElement.h>
#include <Corotational.d/SuperCorotator.h>
#include <Utils.d/dbg_alloca.h>
#include <Driver.d/Mpc.h>

SuperElement::SuperElement(bool _localFlag)
 : eset(0), dsa(0), superCorotator(0), nInternalNodes(0), css(0),
   subElems(0), nSubElems(0), subElemDofs(0), subElemNodes(0),
   nnodes(0), ndofs(0), nn(0), sub_extf(0)
{ 
  localFlag = _localFlag;
}

SuperElement::~SuperElement()
{
  if(eset) delete eset;
  if(dsa) delete dsa;
  if(css) delete css;
  if(subElems) {
    for(int i = 0; i < nSubElems; ++i) delete subElems[i];
    delete [] subElems;
  }
  if(subElemDofs) {
    for(int i = 0; i < nSubElems; ++i) delete [] subElemDofs[i];
    delete [] subElemDofs;
  }
  if(subElemNodes) {
    for(int i = 0; i < nSubElems; ++i) delete [] subElemNodes[i];
    delete [] subElemNodes;
  }
  if(nn) delete [] nn;
  if(sub_extf) {
    for(int i = 0; i < nSubElems; ++i) delete [] sub_extf[i];
    delete [] sub_extf;
  }
}

void
SuperElement::setPreLoad(double load, int &flg)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setPreLoad(load, flg);
}

void
SuperElement::setPressure(double pres, MFTTData *mftt)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setPressure(pres, mftt);
}

void
SuperElement::setFrame(EFrame *frame)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setFrame(frame);
}

void
SuperElement::buildFrame(CoordSet &cs)
{
  if(localFlag) {
    css = new CoordSet(nnodes); // coordinate subset
    for(int i = 0; i < nnodes; ++i) (*css)[i] = cs[nn[i]];
    for(int i = 0; i < nSubElems; ++i) subElems[i]->buildFrame(*css);
  }
  else {
    for(int i = 0; i < nSubElems; ++i) subElems[i]->buildFrame(cs);
  }
}

void 
SuperElement::setProp(StructProp *_prop, bool _myProp) 
{
  if(myProp && prop) delete prop;

  prop = _prop; 
  myProp = _myProp;
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setProp(prop, false);
  makeAllDOFs();
}

void
SuperElement::setCompositeData(int _type, int nlays, double *lData,
                               double *coefs, double *frame)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setCompositeData(_type, nlays, lData, coefs, frame);
}

double *
SuperElement::setCompositeData2(int _type, int nlays, double *lData,
                                double *coefs, CoordSet &cs, double theta)
{
  // compute the material coordinate frame for the first sub-element
  // assuming subElems[0] nodes 0 and 1 are same as superElement nodes 0 and 1
  // which is the case for Compo4NodeShell
  double *frame = subElems[0]->setCompositeData2(_type, nlays, lData, coefs, cs, theta);

  // assign the frame from subElem[0] to all the other subElems
  for(int i = 1; i < nSubElems; ++i) subElems[i]->setCompositeData(_type, nlays, lData, coefs, frame);
  return frame;
}

void
SuperElement::setOffset(double *o)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setOffset(o);
}

void::
SuperElement::setMaterial(NLMaterial *m)
{
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setMaterial(m);
}

int
SuperElement::numInternalNodes()
{
  // nInternalNodes has already been set in makeAllDOFs
  return nInternalNodes;
}

void
SuperElement::setInternalNodes(int *in)
{
  int offset = 0;
  // set sub-element internal nodes
  for(int i = 0; i < nSubElems; ++i) {
    subElems[i]->setInternalNodes(in+offset);
    offset += subElems[i]->numInternalNodes();
  }
  // add super-element internal node numbers to nn
  for(int i = 0; i < offset; ++i) nn[nnodes+i] = in[i];
}

void
SuperElement::renum(int *table)
{
  for(int i = 0; i < numNodes(); ++i) nn[i] = table[nn[i]];
  for(int i = 0; i < nSubElems; ++i) subElems[i]->renum(table);
}

FullSquareMatrix 
SuperElement::stiffness(CoordSet &cs, double *karray, int flg)
{
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  for(int i = 0; i < nSubElems; ++i) {
    double *subkarray = (double *) dbg_alloca(sizeof(double)*subElems[i]->numDofs()*subElems[i]->numDofs());
    FullSquareMatrix subk = subElems[i]->stiffness(cs, subkarray, flg);
    ret.add(subk, subElemDofs[i]);
  }
  return ret;
}

FullSquareMatrix 
SuperElement::massMatrix(CoordSet &cs, double *marray, int cmflg)
{
  FullSquareMatrix ret(numDofs(), marray);
  ret.zero();
  for(int i = 0; i < nSubElems; ++i) {
    double *submarray = (double *) dbg_alloca(sizeof(double)*subElems[i]->numDofs()*subElems[i]->numDofs());
    FullSquareMatrix subm = subElems[i]->massMatrix(cs, submarray, (double)cmflg);
    ret.add(subm, subElemDofs[i]);
  }
  return ret;
}

double 
SuperElement::getMass(CoordSet &cs)
{
  double ret = 0.0;
  for(int i = 0; i < nSubElems; ++i) ret += subElems[i]->getMass(cs);
  return ret;
}

void
SuperElement::getGravityForce(CoordSet &cs, double *gravityAcceleration, Vector &gravityForce,
                              int gravflg, GeomState *geomState)
{
  gravityForce.zero();
  for(int i = 0; i < nSubElems; ++i) {
    Vector subGravityForce(subElems[i]->numDofs());
    subElems[i]->getGravityForce(cs, gravityAcceleration, subGravityForce, gravflg, geomState);
    gravityForce.add(subGravityForce, subElemDofs[i]);
  }
}

void 
SuperElement::getThermalForce(CoordSet &cs, Vector &nodeTemp, Vector &thermalForce,
                              int glflag, GeomState *gs)
{
  thermalForce.zero();
  for(int i = 0; i < nSubElems; ++i) {
    Vector subThermalForce(subElems[i]->numDofs());
    Vector subNodeTemp(nodeTemp, subElems[i]->numNodes(), subElemNodes[i]);
    subElems[i]->getThermalForce(cs, subNodeTemp, subThermalForce, glflag, gs);
    thermalForce.add(subThermalForce, subElemDofs[i]);
  }
}

void 
SuperElement::getIntrnForce(Vector &elForce, CoordSet &cs,
                            double *elDisp, int Index, double *nodeTemp)
{
  elForce.zero();
  Vector nodeTempVec(nodeTemp, numNodes());
  Vector elDispVec(elDisp, numDofs());
  for(int i = 0; i < nSubElems; ++i) {
    Vector subElementForce(subElems[i]->numDofs());

    double *subElementDisp = 0;
    // if available, use the element displacemented from the corotator instead of elDisp (for non-linear)
    if(superCorotator) subElementDisp = superCorotator->getPreviouslyExtractedSubDeformations(i);
    bool delete_flag = false;
    if(!subElementDisp) {
      subElementDisp = new double[subElems[i]->numDofs()]; 
      delete_flag = true;
      for(int j = 0; j < subElems[i]->numDofs(); ++j) subElementDisp[j] = elDisp[subElemDofs[i][j]];
    }

    double *subNodeTemp = 0;
    if(nodeTemp) {
      subNodeTemp = new double[subElems[i]->numNodes()];
      for(int j = 0; j < subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
    }

    subElems[i]->getIntrnForce(subElementForce, cs, subElementDisp, Index, subNodeTemp);

    elForce.add(subElementForce, subElemDofs[i]);

    if(subNodeTemp) delete [] subNodeTemp;
    if(delete_flag) delete [] subElementDisp;
  }
}

void 
SuperElement::getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                          Vector &elDisp, int strInd, int surface,
                          double *nodeTemp, double ylayer, double zlayer, int avgnum)
{
  stress.zero();
  weight.zero();
  for(int i = 0; i < nSubElems; ++i) {
    Vector subElementStress(subElems[i]->numNodes());
    Vector subElementWeight(subElems[i]->numNodes());

    Vector *subElementDisp = 0;
    if(superCorotator) {
      double *subd = superCorotator->getPreviouslyExtractedSubDeformations(i);
      // if available, use the element displacemented from the corotator instead of elDisp (for non-linear)
      if(subd) subElementDisp = new StackVector(subd, subElems[i]->numDofs());
    }
    if(!subElementDisp) subElementDisp = new Vector(elDisp, subElems[i]->numDofs(), subElemDofs[i]);

    double *subNodeTemp = 0;   
    if(nodeTemp) {
      subNodeTemp = new double[subElems[i]->numNodes()];
      for(int j = 0; j < subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
    }

    subElems[i]->getVonMises(subElementStress, subElementWeight, cs, *subElementDisp, strInd,
                             surface, subNodeTemp, ylayer, zlayer, avgnum);

    stress.add(subElementStress, subElemNodes[i]);
    weight.add(subElementWeight, subElemNodes[i]);

    if(subNodeTemp) delete [] subNodeTemp;
    delete subElementDisp;
  }
}

void 
SuperElement::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface, double *nodeTemp)
{
  stress.zero();
  weight.zero();
  for(int i = 0; i < nSubElems; ++i) {
    FullM subElementStress(subElems[i]->numNodes(), 9);
    Vector subElementWeight(subElems[i]->numNodes());

    Vector *subElementDisp = 0;
    if(superCorotator) {
      double *subd = superCorotator->getPreviouslyExtractedSubDeformations(i);
      // if available, use the element displacemented from the corotator instead of elDisp (for non-linear)
      if(subd) subElementDisp = new StackVector(subd, subElems[i]->numDofs());
    }
    if(!subElementDisp) subElementDisp = new Vector(elDisp, subElems[i]->numDofs(), subElemDofs[i]);

    double *subNodeTemp = 0;
    if(nodeTemp) {
      subNodeTemp = new double[subElems[i]->numNodes()];
      for(int j = 0; j < subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
    }

    subElems[i]->getAllStress(subElementStress, subElementWeight, cs, *subElementDisp, strInd,
                              surface, subNodeTemp);

    stress.addrows(subElementStress, subElemNodes[i]);
    weight.add(subElementWeight, subElemNodes[i]);
 
    if(subNodeTemp) delete [] subNodeTemp;
    delete subElementDisp;
  }
}

void 
SuperElement::computeHeatFluxes(Vector &heatflux, CoordSet &cs, Vector &elTemp, int hflInd)
{
  heatflux.zero();
  for(int i = 0; i < nSubElems; ++i) {
    Vector subElementHeatFlux(subElems[i]->numDofs());
    Vector subElementTemp(elTemp, subElems[i]->numNodes(), subElemNodes[i]);
    subElems[i]->computeHeatFluxes(subElementHeatFlux, cs, subElementTemp, hflInd);
    heatflux.add(subElementHeatFlux, subElemDofs[i]);
  }
}

void 
SuperElement::trussHeatFluxes(double &trussflux, CoordSet &cs, Vector &elTemp, int hflInd)
{
  trussflux = 0.0;
  for(int i = 0; i < nSubElems; ++i) {
    double subTrussFlux = 0.0;
    Vector subElementTemp(elTemp, subElems[i]->numNodes(), subElemNodes[i]);
    subElems[i]->trussHeatFluxes(subTrussFlux, cs, subElementTemp, hflInd);
    trussflux += subTrussFlux;
  }
}

void 
SuperElement::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs)
{
  cerr << " *** WARNING: SuperElement::computeDisp(...) is not implemented \n";
}

void
SuperElement::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                        double *resF, GeomState *gs)
{
  cerr << " *** WARNING: SuperElement::getFlLoad(...) is not implemented \n";
}

void
SuperElement::computeTemp(CoordSet &cs, State &state, double gp[2], double *tres)
{
  cerr << " *** WARNING: SuperElement::computeTemp(...) is not implemented \n";
  // need to determine from gauss point values (gp) the appropriate sub element i
  // for which to call subElems[i]->computeTemp(cs, state, subElem_gp, tres)
}

void
SuperElement::getFlFlux(double gp[2], double *flF, double *tresF)
{
  cerr << " *** WARNING: SuperElement::getFlFlux(...) is not implemented \n";
}

void
SuperElement::markDofs(DofSetArray &dsa)
{
  for(int i = 0; i < nSubElems; ++i)
    subElems[i]->markDofs(dsa);
}

int*
SuperElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];
  for(int i = 0; i < nSubElems; ++i) {
    int *subp = new int[subElems[i]->numDofs()];
    subp = subElems[i]->dofs(dsa, subp);
    for(int j = 0; j < subElems[i]->numDofs(); ++j) {
      p[subElemDofs[i][j]] = subp[j];
    }
    if(subp) delete [] subp;
  }
  return p;
}

int 
SuperElement::numDofs()
{
  return ndofs;
}

int 
SuperElement::numNodes()
{
  return nnodes + nInternalNodes;
}

int* 
SuperElement::nodes(int *p)
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < numNodes(); ++i) p[i] = nn[i];
  return p;
}

Corotator*
SuperElement::getCorotator(CoordSet &cs, double *d, int i1, int i2)
{
  superCorotator = new SuperCorotator(this);
  for(int i = 0; i < nSubElems; ++i)
    superCorotator->setSubCorotator(i, subElems[i]->getCorotator(cs, d, i1, i2));
 
  return superCorotator;
}

void 
SuperElement::computePressureForce(CoordSet &cs, Vector &elPressureForce,
                                   GeomState *gs, int cflg)
{
  if(!sub_extf) { // save a copy of the external force for each sub-element
    sub_extf = new double * [nSubElems];
    for(int i=0; i<nSubElems; ++i) sub_extf[i] = new double[subElems[i]->numDofs()];
  }

  elPressureForce.zero();
  for(int i = 0; i < nSubElems; ++i) {
    Vector subElementPressureForce(subElems[i]->numDofs(), sub_extf[i], false);
    subElementPressureForce.zero();
    subElems[i]->computePressureForce(cs, subElementPressureForce, gs, cflg);
    elPressureForce.add(subElementPressureForce, subElemDofs[i]);
  }
/*
  // experimental: adjust for four node shell
  int *dofCount = new int[numDofs()];
  for(i=0; i<numDofs(); ++i) dofCount[i] = 0;
  int totsubnodes = 0;
  for(i=0; i<nSubElems; ++i) {
    totsubnodes += subElems[i]->numNodes();
    for(int j=0; j<subElems[i]->numDofs(); ++j) 
      dofCount[subElemDofs[i][j]] += 1;
  }
  for(i=0; i<numDofs(); ++i) { 
    double temp = elPressureForce[i]*totsubnodes/numNodes()/dofCount[i];
    elPressureForce[i] = temp;
  }
  delete [] dofCount;
*/
}

double *
SuperElement::getMidPoint(CoordSet &cs)
{
  // this is only correct for superelements in which all subelements have equal area
  double *midPoint = new double[3];
  midPoint[0] = midPoint[1] = midPoint[2] = 0.0;
  for(int i = 0; i < nSubElems; ++i) {
    double *subElemMidPoint = subElems[i]->getMidPoint(cs);
    midPoint[0] += subElemMidPoint[0];
    midPoint[1] += subElemMidPoint[1];
    midPoint[2] += subElemMidPoint[2];
    delete [] subElemMidPoint;
  }
  midPoint[0] /= nSubElems;
  midPoint[1] /= nSubElems;
  midPoint[2] /= nSubElems;
  return midPoint;
}

double * 
SuperElement::getCompositeData(int nl)
{
  return subElems[0]->getCompositeData(nl); // all sub-elements should have same composite data
}

double *
SuperElement::getCompositeFrame()
{
  return subElems[0]->getCompositeFrame(); // all sub-elements should have same composite frame
}

int
SuperElement::getCompositeLayer()
{
  return subElems[0]->getCompositeLayer(); // all sub-elements should have same composite layer
}

int
SuperElement::dim()
{
  int ret = 0;
  for(int i = 0; i < nSubElems; ++i) 
    if(subElems[i]->dim() > ret) ret = subElems[i]->dim();
  return ret;
}

void 
SuperElement::addFaces(PolygonSet *pset)
{
  cerr << " *** WARNING: SuperElement::addFaces(...) is not implemented \n";
}

bool 
SuperElement::isSafe()
{
  for(int i = 0; i < nSubElems; ++i) 
    if(!subElems[i]->isSafe()) return false;
  return true;
}

bool 
SuperElement::isRotMidSideNode(int iNode)
{
  for(int i = 0; i < nSubElems; ++i) {
    for(int j = 0; j < subElems[i]->numNodes(); ++j) {
      if(subElemNodes[i][j] == iNode) {
         if(subElems[i]->isRotMidSideNode(j)) return true; 
      }
    }
  }
  return false;
}

bool 
SuperElement::isMpcElement()
{
  // return true if one of the sub elements is a mpc element
  for(int i = 0; i < nSubElems; ++i)
    if(subElems[i]->isMpcElement()) return true;
  return false;
}

/*
bool
SuperElement::isRigidMpcElement(const DofSet &dset, bool forAllNodes)
{
  if(forAllNodes) {
	  // return true if one of the sub elements is a rigid mpc element
	  for(int i=0; i<nSubElems; ++i)
	    if(!subElems[i]->isRigidMpcElement(dset, forAllNodes))
	      return false;
	  return true;
  }
  // return true if one of the sub elements is a rigid mpc element
  for(int i=0; i<nSubElems; ++i)
    if(subElems[i]->isRigidMpcElement(dset, forAllNodes)) return true;
  return false;
}
*/

bool
SuperElement::isConstraintElement()
{
  // return true if one of the sub elements is a rigid mpc element
  for(int i = 0; i < nSubElems; ++i)
    if(subElems[i]->isConstraintElement()) return true;
  return false;
}

int
SuperElement::getMassType()
{
  return subElems[0]->getMassType();
}

int
SuperElement::getNumMPCs()
{
  int ret = 0;
  for(int i = 0; i < nSubElems; ++i) ret += subElems[i]->getNumMPCs();
  return ret;
}

LMPCons**
SuperElement::getMPCs()
{
  LMPCons** ret = new LMPCons * [getNumMPCs()];
  int k = 0;
  for(int i = 0; i < nSubElems; ++i) {
    LMPCons** submpcs = subElems[i]->getMPCs();
    for(int j = 0; j < subElems[i]->getNumMPCs(); ++j) {
      ret[k] = new LMPCons(*submpcs[j]);
      k++;
    }
    delete [] submpcs;
  }
  return ret;
}

void
SuperElement::setGlNum(int gn, int sn)
{ 
  glNum = gn;
  subNum = sn;
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setGlNum(glNum, i);
}

int
SuperElement::numStates()
{
  int ns = 0;
  for(int i = 0; i < nSubElems; ++i) ns += subElems[i]->numStates();
  return ns;
}

void
SuperElement::makeAllDOFs()
{
  nInternalNodes = 0;
  if(localFlag) {
    // count and locally number the internal nodes
    for(int i = 0; i < nSubElems; ++i) {
      int k = subElems[i]->numInternalNodes();
      if(k > 0) {
        int *in = new int[k];
        for(int j = 0; j < k; ++j) in[j] = nnodes + (nInternalNodes++);
        subElems[i]->setInternalNodes(in);
        delete [] in;
      }
    }
    // resize nn to allow for internal nodes (if any)
    if(nInternalNodes > 0) {
      int *new_nn = new int[nnodes+nInternalNodes];
      for(int i = 0; i < nnodes; ++i) new_nn[i] = nn[i];
      for(int i = 0; i < nInternalNodes; ++i) new_nn[nnodes+i] = -1;
      if(nn) delete [] nn;
      nn = new_nn;
    }

    // make the sub-to-super node numbering maps
    subElemNodes = new int * [nSubElems];
    for(int i = 0; i < nSubElems; ++i) {
      subElemNodes[i] = new int[subElems[i]->numNodes()];
      subElems[i]->nodes(subElemNodes[i]);
    }

    // make the element set
    eset = new Elemset(nSubElems);
    eset->setMyData(false);
    for(int i = 0; i < nSubElems; ++i) eset->elemadd(i, subElems[i]);

    // make the sub-to-super dof numbering maps
    dsa = new DofSetArray(nnodes+nInternalNodes, *eset);
    subElemDofs = new int * [nSubElems];
    for(int i = 0; i < nSubElems; ++i) {
      subElemDofs[i] = new int[subElems[i]->numDofs()];
      subElems[i]->dofs(*dsa, subElemDofs[i]);
    }
    ndofs = dsa->size();

    // propogate the node numbering of this element down to the sub elements
    for(int i = 0; i < nSubElems; ++i) subElems[i]->renum(nn);
    localFlag = false;
  }
  else {
    for(int i = 0; i < nSubElems; ++i) nInternalNodes += subElems[i]->numInternalNodes();
  }
}
