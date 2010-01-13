#include <Element.d/SuperElement.h>
#include <Corotational.d/SuperCorotator.h>
#include <Utils.d/dbg_alloca.h>

// PJSA: note to self, check all the vector/matrix adds - perhaps some should be averages or overwrites

void
SuperElement::renum(int *table)
{
  int i;
  for(i=0; i<nnodes; ++i) nn[i] = table[nn[i]];
  for(i=0; i<nSubElems; ++i) subElems[i]->renum(table);
}

void
SuperElement::buildCorotator(CoordSet &cs) 
{
  if(!superCorotator) {
    superCorotator = new SuperCorotator(this);
    int i;
    for(i=0; i<nSubElems; ++i) 
      subElems[i]->buildCorotator(cs);
  }
}

void 
SuperElement::setProp(StructProp *p, bool _myProp) 
{
  if(myProp && prop) {
    delete prop;
    prop = 0;
  }

  prop = p; 
  myProp = _myProp;
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->setProp(p, false);
}

void 
SuperElement::setPreLoad(double load, int &flg)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->setPreLoad(load, flg);
}

void
SuperElement::setFrame(EFrame *frame)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->setFrame(frame);
}

void
SuperElement::buildFrame(CoordSet& cs)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->buildFrame(cs);
}

void 
SuperElement::setOffset(double *o)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->setOffset(o);
}

void 
SuperElement::setCompositeData(int _type, int nlays, double *lData,
                               double *coefs, double *frame)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->setCompositeData(_type, nlays, lData, coefs, frame);
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
  int i;
  for(i=1; i<nSubElems; ++i) subElems[i]->setCompositeData(_type, nlays, lData, coefs, frame);
  return frame;
}


FullSquareMatrix 
SuperElement::stiffness(CoordSet &cs, double *karray, int flg)
{
  FullSquareMatrix ret(numDofs(), karray);
  ret.zero();
  int i;
  for(i=0; i<nSubElems; ++i) {
    int size = sizeof(double)*subElems[i]->numDofs()*subElems[i]->numDofs();
    double *subkarray = (double *) dbg_alloca(size);
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
  int i;
  for(i=0; i<nSubElems; ++i) {
    int size = sizeof(double)*subElems[i]->numDofs()*subElems[i]->numDofs();
    double *submarray = (double *) dbg_alloca(size);
    FullSquareMatrix subm = subElems[i]->massMatrix(cs, submarray, (double)cmflg);
    ret.add(subm, subElemDofs[i]);
  }
  return ret;
}

double 
SuperElement::getMass(CoordSet &cs)
{
  double ret = 0.0;
  int i;
  for(i=0; i<nSubElems; ++i) ret += subElems[i]->getMass(cs);
  return ret;
}

void
SuperElement::getGravityForce(CoordSet &cs, double *gravityAcceleration, Vector &gravityForce,
                              int gravflg, GeomState *geomState)
{
  gravityForce.zero();
  int i;
  for(i=0; i<nSubElems; ++i) {
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
  int i;
  for(i=0; i<nSubElems; ++i) {
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
  int i, j;
  Vector nodeTempVec(nodeTemp, numNodes());
  Vector elDispVec(elDisp, numDofs());
  for(i=0; i<nSubElems; ++i) {
    Vector subElementForce(subElems[i]->numDofs());

    double *subElementDisp = 0;
    // if available, use the element displacemented from the corotator instead of elDisp (for non-linear)
    if(superCorotator) subElementDisp = superCorotator->getPreviouslyExtractedSubDeformations(i);
    bool delete_flag = false;
    if(!subElementDisp) {
      subElementDisp = new double[subElems[i]->numDofs()]; 
      delete_flag = true;
      for(j=0; j<subElems[i]->numDofs(); ++j) subElementDisp[j] = elDisp[subElemDofs[i][j]];
    }

    double *subNodeTemp = 0;
    if(nodeTemp) {
      subNodeTemp = new double[subElems[i]->numNodes()];
      for(j=0; j<subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
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
  int i,j;
  for(i=0; i<nSubElems; ++i) {
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
      for(j=0; j<subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
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
  stress.zero(); // this is size nnodes * 9;
  weight.zero();
  int i, j;
  for(i=0; i<nSubElems; ++i) {
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
      for(j=0; j<subElems[i]->numNodes(); ++j) subNodeTemp[j] = nodeTemp[subElemNodes[i][j]];
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
  int i;
  for(i=0; i<nSubElems; ++i) {
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
  int i;
  for(i=0; i<nSubElems; ++i) {
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
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->markDofs(dsa);
}

int*
SuperElement::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];
  int i, j;
  for(i=0; i<nSubElems; ++i) {
    int *subp = new int[subElems[i]->numDofs()];
    subp = subElems[i]->dofs(dsa, subp);
    for(j=0; j<subElems[i]->numDofs(); ++j) {
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
  return nnodes;
}

int* 
SuperElement::nodes(int *p)
{
  if(p == 0) p = new int[numNodes()];
  int i;
  for(i=0; i<numNodes(); ++i) p[i] = nn[i];
  return p;
}

Corotator*
SuperElement::getCorotator(CoordSet &cs, double *d, int i1, int i2)
{
  if(!superCorotator) superCorotator = new SuperCorotator(this);
  int i;
  for(i=0; i<nSubElems; ++i) {
    superCorotator->setSubCorotator(i, subElems[i]->getCorotator(cs, d, i1, i2));
  }
  return superCorotator;
}

void 
SuperElement::setPressure(double pres)
{
  for(int i=0; i<nSubElems; ++i) {
    subElems[i]->setPressure(pres);
  }
}

void 
SuperElement::computePressureForce(CoordSet &cs, Vector &elPressureForce,
                                   GeomState *gs, int cflg)
{
  elPressureForce.zero();
  int i;
  for(i=0; i<nSubElems; ++i) {
    Vector subElementPressureForce(subElems[i]->numDofs());
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
  int i;
  for(i=0; i<nSubElems; ++i) {
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

double 
SuperElement::getMoment(Vector& force, CoordSet& cs, int node, int idir)
{
  cerr << " *** WARNING: SuperElement::getMoment(...) is not implemented \n";
  return 0.0;
}

int
SuperElement::dim()
{
  int ret = 0;
  int i;
  for(i=0; i<nSubElems; ++i) 
    if(subElems[i]->dim() > ret) ret = subElems[i]->dim();
  return ret;
}

void 
SuperElement::addFaces(PolygonSet *pset)
{
  cerr << " *** WARNING: SuperElement::addFaces(...) is not implemented \n";
}

int 
SuperElement::numInternalNodes()
{
  int ret = 0;
  int i;
  for(i=0; i<nSubElems; ++i) ret += subElems[i]->numInternalNodes();
  return ret;
}

void 
SuperElement::setInternalNodes(int *in)
{
  int offset = 0;
  int i;
  // set sub-element internal nodes
  for(i=0; i<nSubElems; ++i) {
    subElems[i]->setInternalNodes(in+offset);
    offset += subElems[i]->numInternalNodes();
  }
  // set super-element internal nodes
  for(i=0; i<offset; ++i) nn[nnodes-offset+i] = in[i];
}

bool 
SuperElement::isSafe()
{
  int i;
  for(i=0; i<nSubElems; ++i) 
    if(!subElems[i]->isSafe()) return false;
  return true;
}

bool 
SuperElement::isRotMidSideNode(int iNode)
{
  int i,j;
  for(i=0; i<nSubElems; ++i) {
    for(j=0; j<subElems[i]->numNodes(); ++j) {
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
  int i;
  for(i=0; i<nSubElems; ++i)
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
  for(int i=0; i<nSubElems; ++i)
    if(subElems[i]->isConstraintElement()) return true;
}

/*
void 
SuperElement::computeMPCs(CoordSet &cs)
{
  int i;
  for(i=0; i<nSubElems; ++i) subElems[i]->computeMPCs(cs);
}
*/

SuperElement::~SuperElement()
{
/*
  int i;
  for(i=0; i<nSubElems; ++i) { 
    if(subElems[i]) delete subElems[i];
    if(subElemNodes[i]) delete [] subElemNodes;
   if(subElemDofs[i]) delete [] subElemDofs;
  }
  if(subElems) delete [] subElems;
  if(subElemNodes) delete [] subElemNodes;
  if(subElemDofs) delete [] subElemDofs;
  if(nn) delete [] nn;
*/
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
  int i;
  for(i=0; i<nSubElems; ++i) ret += subElems[i]->getNumMPCs();
  return ret;
}

LMPCons**
SuperElement::getMPCs()
{
  LMPCons** ret = new LMPCons * [getNumMPCs()];
  int i,j,k=0;
  for(i=0; i<nSubElems; ++i) {
    LMPCons** submpcs = subElems[i]->getMPCs();
    for(j=0; j<subElems[i]->getNumMPCs(); ++j) ret[k++] = submpcs[j];
  }
  return ret;
}

void
SuperElement::initialize(int l, int* _nn)
{
  // this function is designed to be called in the constructor of a super element
  // after the sub elements have been instantiated. See for example Element.d/Joint.d/RigidJoint.C
  // l is the number of nodes, excluding internal nodes
  // _nn is an array of dimension l containing the node numbers

  // make the element set
  Elemset eset; eset.setMyData(false);
  for(int i = 0; i < nSubElems; ++i) eset.elemadd(i, subElems[i]);

  // count and number the internal nodes
  int m = 0;
  for(int i = 0; i < nSubElems; ++i) {
    int k = subElems[i]->numInternalNodes();
    int *in = new int[k];
    for(int j = 0; j < k; ++j) in[j] = l + (m++);
    subElems[i]->setInternalNodes(in);
    delete [] in;
  }
  nnodes = l + m;

  // get the sub-to-super node numbering maps
  subElemNodes = new int * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) { 
    subElemNodes[i] = new int[subElems[i]->numNodes()];
    subElems[i]->nodes(subElemNodes[i]); 
  }
  //for(int i=0; i<nSubElems; ++i) { cerr << "subElemNodes[" << i << "] = "; for(int j=0; j<subElems[i]->numNodes(); ++j) cerr << subElemNodes[i][j] << " "; cerr << endl; }

  // number the sub-to-super dof numbering maps
  DofSetArray dsa(nnodes, eset);
  subElemDofs = new int * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) { 
    subElemDofs[i] = new int[subElems[i]->numDofs()]; 
    subElems[i]->dofs(dsa, subElemDofs[i]); 
  }
  ndofs = dsa.size();

  // renumber sub element nodes to global node numbering
  int *new_nn = new int[nnodes];
  for(int i = 0; i < l; ++i) new_nn[i] = _nn[i]; for(int i = l; i < nnodes; ++i) new_nn[i] = -1;
  for(int i = 0; i < nSubElems; ++i) subElems[i]->renum(new_nn);
  if(nn) delete [] nn; nn = new_nn;

  for(int i = 0; i < nSubElems; ++i) subElems[i]->setGlNum(-1);  
}
