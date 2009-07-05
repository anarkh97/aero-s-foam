#include <stdio.h>
#include <stdlib.h>
#include <sstream>    //CRW
#include <Utils.d/BlockAlloc.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/Mpc.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/Connectivity.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Utils.d/ModeData.h>
#include <Dec.d/Decomp.d/Decomp.h>
#include <Driver.d/MultiFront.h>
#include <Driver.d/Header.h>
#ifndef WINDOWS
#include <dlfcn.h>
#endif
#include <map>
#include <utility>
using std::map;
#include <list>
#include <Sfem.d/Sfem.h>

EFrameData null_eframe;
StructProp null_sprops;
Attrib null_attrib;
BCond null_bcond;
CoefData null_coef;
OutputInfo emptyInfo;
extern ModeData modeData;
#ifndef SALINAS
extern Sfem *sfem;
#endif
//----------------------------------------------------------------------

GeoSource::GeoSource(int iniSize) : oinfo(emptyInfo, iniSize), nodes(iniSize*16), elemSet(iniSize*16),
   layInfo(0, iniSize), coefData(0, iniSize), layMat(0, iniSize), efd(null_eframe, iniSize), cframes(0, iniSize)
{
  decJustCalled=false;
  exitAfterDec=false;
  nGlobNodes = 0;
  numOutInfo = 0;
  outLimit = -1;
  numNodalOutput = 0;
  outputNodes = 0;
  outNodeIndex = 0;
  headLen = 0;
  curCluster = -1;

  nElem = 0;
  nElemFluid = 0;
  phantomFlag = 0;  // 0 => no phantom elems
  numClusters = 0;

  elemTypeNumNodesMap = 0;
  na = 0;
  namax = 0;
  numEframes = 0;
  numCframes = 0;
  numLayInfo = 0;
  numCoefData = 0;
  numLayMat = 0;
  prsflg = 0;
  prlflg = 0;

  constpflg = 0;
  constqflg = 0;
  maxGlobNode = 0;
  maxClusNode = 0;
  numProps = 0;

  // set file names to NULL
  mapName = (char *) "CPUMAP";  // default cpu map
  matchName = NULL;

  // initialize bc's
  numTextDirichlet = 0;
  numTextNeuman = 0;
  numDirichlet = 0;
  numDirichletFluid = 0; //ADDED FOR HEV PROBLEM, EC, 20070820
  numNeuman = 0;
  numConvBC = 0;
  numRadBC = 0;
  numIDis = 0;
  numIDis6 = 0;
  numIVel = 0;
  numITemp = 0;
  numDampedModes = 0;
  numComplexDirichlet = 0;
  numComplexNeuman = 0;
  numSurfaceDirichlet = 0;
  numSurfaceNeuman = 0;
  numSurfacePressure = 0;

  // PITA
  // Initial seed conditions
  PitaIDis6 = PitaIVel6 = 0;
  numPitaIDis6   = 0;
  numTSPitaIDis6 = 0;
  numPitaIVel6   = 0;
  numTSPitaIVel6 = 0;
  newStep0 = false;

  subToCPU = 0;
  cpuToSub = 0;
  cpuToCPU = 0;
  subToNode = 0;
  subToElem = 0;
  subToSub = 0;

  // init match data
  numMatchData = 0;
  matchData = 0;
  numGapVecs = 0;
  gapVec = 0;

  cinfo = new ControlInfo;
  claw = 0;

  isShift = false;
  shiftV = 0.0;

  dbc = 0;
  dbcFluid = 0;
  nbc = 0;
  textDBC = dbc = textNBC = nbc = cvbc = iDis = iDis6 = iVel = iTemp = modalDamping = 0;
  cdbc = cnbc = 0;
  surface_dbc = surface_nbc = surface_pres = 0;

  maxattrib = -1; // PJSA
  optDec = 0;

  numInternalNodes = 0;
  allNumClusElems = 0;
  binaryInput = false;
  binaryInputControlLeft = false;
#ifdef DISTRIBUTED
  binaryOutput = true;
#else
  binaryOutput = false;
#endif
  subToClus = 0;

  mpcDirect = false;

  initMRatio();
}

//----------------------------------------------------------------------

void GeoSource::cleanUp()
{
  nodes.deleteNodes();
  elemSet.deleteElems();
  nElem = 0;
  nElemFluid = 0; //ADDED FOR HEV PROBLEM, EC, 20070820
}

//----------------------------------------------------------------------

GeoSource::~GeoSource()
{
  /* not fully implemented */
  if(cpuToSub) delete cpuToSub;
  if(cinfo) delete cinfo;
  if(outputNodes) delete [] outputNodes;
  if(outNodeIndex) delete [] outNodeIndex;
  if(headLen) delete [] headLen;
  if(subToClus) delete subToClus;
  // claw is deleted by domain
}

//----------------------------------------------------------------------

int GeoSource::addNode(int nd, double xyz[3])
{
  nodes.nodeadd(nd, xyz);
  nGlobNodes++;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addElem(int en, int type, int nn, int *nodeNumbers)
{
#ifndef SALINAS
  elemSet.elemadd(en, type, nn, nodeNumbers);
#else
  cerr << "*** ERROR: GeoSource::addElem(...) not included in Salinas library \n";
#endif
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addMat(int nmat, StructProp &p)
{
  if (numProps < nmat+1) // attempt to get numProps right -- Julien & Thomas
    numProps = nmat+1;

  p.soundSpeed = omega()/complex<double>(p.kappaHelm, p.kappaHelmImag); // PJSA 1-15-08

  sProps[nmat] = p;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addLay(int l, LayInfo *linfo)
{
  if (numLayInfo <= l)
    numLayInfo = l+1;

  layInfo[l] = linfo;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addCoefInfo(int cin, CoefData &cdata)  {

  if (numCoefData <= cin)
    numCoefData = cin+1;

  CoefData *cd  = new CoefData(cdata);
  coefData[cin] = cd;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addLayMat(int m, double *d)
{
  if (numLayMat <= m)
    numLayMat = m+1;

  LayMat *lm = new LayMat(m, d);
  layMat[m] = lm;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::setAttrib(int n, int a, int ca, int cfrm, double ctheta)
{
  if(a>maxattrib) maxattrib = a; // PJSA
  attrib[na].nele = n;
  attrib[na].attr = a;
  attrib[na].cmp_attr = ca;
  attrib[na].cmp_frm  = cfrm;
  attrib[na].cmp_theta = ctheta;
  na++;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::setFrame(int el, double *data)
{
  efd[numEframes].elnum = el;
  int i,j;
  for(i = 0; i < 3; ++i)
   for(j = 0; j < 3; ++j)
     efd[numEframes].frame[i][j] = data[i*3+j];

  numEframes++;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addCFrame(int fn, double *f)  {

  if (fn >= numCframes)
    numCframes = fn+1;

  cframes[fn] = new double[9];

  for (int i = 0; i < 9; ++i)
    cframes[fn][i] = f[i];

  return 0;
}

//----------------------------------------------------------------------

void GeoSource::addMpcElements(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
  //std::cerr << "MLXX Adding " << numLMPC << " MPC Elements" << std::endl;
#ifndef SALINAS
  // PJSA: convert rigid elements into mpcs -- used for GRBM
  // 1. find largest lmpcnum
  int lmpcnum = 0;
  bool print_flag = true;
  for(int i=0; i < numLMPC; ++i)
    if(lmpc[i]->lmpcnum > lmpcnum) lmpcnum = lmpc[i]->lmpcnum;
  // 2. loop through all elements & if element is a rigid element then compute lmpc & add to global list
  int nEle = elemSet.last();
  for(int i = 0; i < na; ++i) {
    Element *ele = elemSet[ attrib[i].nele ];
    if((ele != 0) && (ele->isRigidMpcElement())) {
      if(print_flag) { fprintf(stderr," ... Converting RigidMpcElements to LMPCs \n"); print_flag = false; }
      ele->setProp(&sProps[attrib[i].attr]);  // PJSA 9-18-06 rigid springs need prop
      ele->computeMPCs(nodes,lmpcnum);
    }
  }
  numLMPC = domain->getNumLMPC(); // update to include RigidMpcElements

  if(numLMPC) {

    fprintf(stderr," ... Converting LMPCs to MpcElements \n");
    for(int i = 0; i < numLMPC; ++i) {
      elemSet.mpcelemadd(nEle, lmpc[i]);
      nEle++;
    }
    // XXXX still needed for eigen GRBM lmpc.deleteArray(); domain->setNumLMPC(0);
  }
#else
  cerr << "*** ERROR: GeoSource::addMpcElements(...) not included in Salinas library \n";
#endif
}

/** Order the terms in MPCs so that the first term can be directly written in terms of the others */
void GeoSource::makeDirectMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc) {
  std::cerr << "Making direct MPCs" << std::endl;
  int lmpcnum = 0;
  bool print_flag = true;
  for(int i=0; i < numLMPC; ++i)
    if(lmpc[i]->lmpcnum > lmpcnum) lmpcnum = lmpc[i]->lmpcnum;
  // 2. loop through all elements & if element is a rigid element then compute lmpc & add to global list
  int nEle = elemSet.last();
  Elemset rigidSet;
  int nRigid = 0;
  for(int i = 0; i < na; ++i) {
    Element *ele = elemSet[ attrib[i].nele ];
    if((ele != 0) && (ele->isRigidMpcElement())) {
      if(print_flag) { fprintf(stderr," ... Converting RigidMpcElements to LMPCs \n"); print_flag = false; }
      ele->setProp(&sProps[attrib[i].attr]);  // PJSA 9-18-06 rigid springs need prop
 // Delay this 'til later     ele->computeMPCs(nodes,lmpcnum);
    }
  }
  int n101 = 0;
  for(int i = 0; i < nEle; ++i) {
	 Element *ele = elemSet[i];
	 if(ele != 0 && ele->isRigidMpcElement())
         rigidSet.elemadd(nRigid++,ele);
  }
  std::cerr << "Number of rigid elements: " << nRigid << " versus "<< nEle << std::endl;
  rigidSet.collapseRigid6();
  nRigid = rigidSet.last();
  for(int i = 0; i < nRigid; ++i) {
	  rigidSet[i]->computeMPCs(nodes,lmpcnum);
  }
  numLMPC = domain->getNumLMPC(); // update to include RigidMpcElements
  fprintf(stderr," New number of MPCs: %d\n", numLMPC);
  if(numLMPC) {
    std::map<pair<int,int>, int> tCount;
    for(int i=0; i < numLMPC; ++i) {
      for(int j = 0; j < lmpc[i]->nterms; ++j) {
        if(lmpc[i]->terms[j].isNull())
          continue;
        pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
        std::map<pair<int,int>, int>::iterator it =
          tCount.find(p);
        if(it == tCount.end())
          tCount[p] = 1;
        else
          it->second++;
      }
    }

    for(int i=0; i < numLMPC; ++i) {
      int j;
      for(j = 0; j < lmpc[i]->nterms; ++j) {
        if(lmpc[i]->terms[j].isNull())
          continue;
        pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
        if(tCount[p] == 1) {
          if(j != 0) {
            LMPCTerm lmTerm = lmpc[i]->terms[j];
            lmpc[i]->terms[j] = lmpc[i]->terms[0];
            lmpc[i]->terms[0] = lmTerm;
          }
          break;
        }
      }
      if(j == lmpc[i]->nterms) {
        fprintf(stderr, "ERROR: MPC %d could not be written independently. Please re-organize the MPCs!\n", i+1);
        for(j = 0; j < lmpc[i]->nterms; ++j) {
            pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
        	fprintf(stderr, "%d %d %e - %d, ",
        			lmpc[i]->terms[j].nnum+1, lmpc[i]->terms[j].dofnum+1,
        			lmpc[i]->terms[j].coef.r_value, tCount[p]);
        }
      }
    }

  }
  std::cerr<< "NElement = " << nEle << std::endl;
}

//----------------------------------------------------------------------

void GeoSource::addFsiElements(int numFSI, ResizeArray<LMPCons *> &fsi)
{
#ifndef SALINAS
  if(numFSI) filePrint(stderr," ... Adding FSI Elements \n");
  int nEle = elemSet.last();
  int i, j;
// JLchange:
//  for(i = 0; i < numFSI; ++i) {
//    elemSet.fsielemadd(nEle, fsi[i]);
//    nEle++;
//  }
// Instead of adding one fsi to ONE element, add each term (between one fluid node and one
// structure node) to a single element.
  for(i = 0; i < numFSI; ++i) {
    int fluidNode = fsi[i]->lmpcnum;
    for(j=0; j< (fsi[i])->nterms; j++) {
      LMPCons *thisFsi = new LMPCons(fluidNode, 0.0);
      LMPCTerm thisLmpcTerm((fsi[i])->terms[j], 1.0);
      thisFsi->addterm(&(thisLmpcTerm));
      elemSet.fsielemadd(nEle, thisFsi);
      nEle++;
    }
  }

//  fsi.deleteArray(); domain->setNumFSI(0); // DEBUG fsi_element
#else
  cerr << "*** ERROR: GeoSource::addFsiElements(...) not included in Salinas library \n";
#endif
}

//----------------------------------------------------------------------
void GeoSource::duplicateFilesForPita(int localNumSlices, const int* sliceRankSet)
{
  int initialNumFiles = numOutInfo;
  std::vector<char*> initialFileNameSet(initialNumFiles);
  OutputInfo oI;

  if (localNumSlices > 0)
  {
    for (int j = 0; j < initialNumFiles; ++j)
    {
      oinfo[j].timeSliceRank = sliceRankSet[0];
      initialFileNameSet[j] = oinfo[j].filename;
    }
    timeSliceOutputFiles.insert(std::make_pair(sliceRankSet[0], std::pair<int, int>(0, initialNumFiles)));
  }

  // If necessary, create a whole new set of output files
  // for each time-slice beyond the first on the local CPU
  // and initialize the timeSliceRank tags
  for (int i = 1; i < localNumSlices; ++i)
  {
    int firstRequest = numOutInfo;
    for (int j = 0; j < initialNumFiles; ++j)
    {
      oI.copyParam(oinfo[j]);
      oI.timeSliceRank = sliceRankSet[i];
      addOutput(oI);
    }
    int lastRequest = numOutInfo;
    timeSliceOutputFiles.insert(std::make_pair(sliceRankSet[i], make_pair(firstRequest, lastRequest)));
  }

  // Append the time-slice rank to the output file names
  for (int i = 0; i < numOutInfo; ++i)
  {
    std::stringstream s;
    s << oinfo[i].filename << '.' << oinfo[i].timeSliceRank;
    const std::string & newFileNameString = s.str();

    int newFileNameLength = newFileNameString.size();
    char* newFileName = (char*) malloc(sizeof(char) * (newFileNameLength + 1));
    newFileNameString.copy(newFileName, newFileNameLength);
    newFileName[newFileNameLength] = '\0';

    oinfo[i].filename = newFileName;
  }

  for (int i = 0; i < initialNumFiles; ++i)
  {
    if (initialFileNameSet[i] != 0)
    {
      free(initialFileNameSet[i]);
      initialFileNameSet[i] = 0;
    }
  }
}


//----------------------------------------------------------------------
void GeoSource::setUpData()
{
  int iElem;

  // Setup the internal nodes
  int lastNode = numNodes = nodes.last();
  int nMaxEle = elemSet.size();
  for(iElem = 0; iElem < nMaxEle; ++iElem)
    {
      Element *ele = elemSet[iElem];
      if(ele == 0) continue;
      int nIN = ele->numInternalNodes();
      int nn[125];
      if(nIN > 125) {
	fprintf(stderr, " *** ERROR: Code was not developed for elements with more "
		"than 125 internal nodes");
	exit(-1);
      }
      for(int i = 0; i < nIN; ++i)
	nn[i] = lastNode++;
      ele->setInternalNodes(nn);
      localToGlobalElementsNumber.push_back(iElem);
  }
  numInternalNodes = lastNode-numNodes;

  // Set up element frames
  for (int iFrame = 0; iFrame < numEframes; iFrame++)  {
    Element *ele = elemSet[efd[iFrame].elnum];
    if(ele == 0) {
      // fprintf(stderr," *** WARNING: Frame was found for non-existent"
      //                " element %d \n", efd[iFrame].elnum+1);
    }
    else ele->setFrame(&(efd[iFrame].frame));
  }

  // Set up element attributes
  int i;

  SolverInfo sinfo = domain->solInfo();
  if((na == 0) && (sinfo.probType != SolverInfo::Top) && (sinfo.probType != SolverInfo::Decomp)) { // PJSA
    fprintf(stderr," **************************************\n");
    fprintf(stderr," *** ERROR: ATTRIBUTES not defined  ***\n");
    fprintf(stderr," **************************************\n");
    exit(-1);
  }
  //if((na != nMaxEle) && ((sinfo.probType == SolverInfo::Top) || (sinfo.probType == SolverInfo::Decomp))) { // PJSA
  if(na != nMaxEle) {
    // check for elements with no attribute, and add dummy properties
    bool *hasAttr = (bool *)dbg_alloca(nMaxEle*sizeof(bool));
    for(i = 0; i < nMaxEle; ++i) hasAttr[i] = false;
    for(i = 0; i < na; ++i) hasAttr[attrib[i].nele] = true;

    int dattr = maxattrib+1;  // dummy attribute
    bool topdec = (sinfo.probType == SolverInfo::Top || sinfo.probType == SolverInfo::Decomp); // top & decomp should work without attributes defined
    for(i = 0; i < nMaxEle; ++i) {
      if(elemSet[i] && !hasAttr[i]) {
        if(topdec || elemSet[i]->isConstraintElement()) setAttrib(i,dattr);
        else cerr << " *** WARNING: Element " << i+1 << " has no attribute defined\n";
      }
    }
  }

  // create map for used material properties
  // & compute average E and nu (for coupled_dph)
//  std::map<int, StructProp *> usedMat;
  int structure_element_count = 0;
  int fluid_element_count = 0;
  global_average_E = 0.0;
  global_average_nu = 0.0;
  global_average_rhof = 0.0;
  for(i = 0; i < na; ++i) {
    Element *ele = elemSet[ attrib[i].nele ];

    // Check if element exists
    if (ele == 0) {
       fprintf(stderr," *** WARNING: Attribute was found for"
                      " non existent element %d \n",attrib[i].nele+1);
      continue;
    }

    if(attrib[i].attr < -1) { // phantom elements
      phantomFlag = 1;
      ele->setProp(0);
    }
    else if(elemSet[attrib[i].nele]->isConstraintElement()) { // rigid/mpc/fsi
      ele->setProp(new StructProp(), true);
    }
    else  {
/* PJSA 1-16-08 this doesn't make sense to me
      StructProp *prop;
      std::map<int, StructProp *>::iterator it = usedMat.find(attrib[i].attr);
      if(it == usedMat.end()) {
        if((!sProps[attrib[i].attr].isReal) && (sinfo.probType != SolverInfo::Top)
            && !domain->solInfo().isAcousticHelm() && (sinfo.probType != SolverInfo::Decomp)
            && !elemSet[attrib[i].nele]->isConstraintElement()) {
          fprintf(stderr, " *** WARNING: The material for element %d does not exist\n",
                  attrib[i].nele+1);
        }
        // prop = new StructProp(sProps[attrib[i].attr]);
        prop = &sProps[attrib[i].attr]; // PJSA: for helmsweep, helmCoef is updated in sProps[0] so don't copy!!
                                        // check with Thuan why he decided to make a copy of the material here
        usedMat[attrib[i].attr] = prop;
      }
      else
        prop = it->second;
*/
      SPropContainer::iterator it = sProps.find(attrib[i].attr);
      if(it == sProps.end()) { // in this case the material does not exist (this is only permitted for -t and --dec)
        if((sinfo.probType != SolverInfo::Top) && (sinfo.probType != SolverInfo::Decomp))
          fprintf(stderr, " *** ERROR: The material for element %d does not exist\n",
                  attrib[i].nele+1);
      }
      else {
        StructProp *prop = &(it->second);

        // compute global average structural and fluid properties
        if(! dynamic_cast<HelmElement *>(elemSet[attrib[i].nele])) { // not a fluid element
          global_average_E += prop->E;
          global_average_nu += prop->nu;
          structure_element_count++;
        } else {
          global_average_rhof += prop->rho;
          fluid_element_count++;
        }

        ele->setProp(prop);
      }
    }

    ele->buildFrame(nodes);

    if(attrib[i].cmp_attr >= 0) {
      if(coefData[attrib[i].cmp_attr] != 0) {
        if(attrib[i].cmp_frm > -1) { // user input cframe
          ele->setCompositeData(1, 0, 0, coefData[attrib[i].cmp_attr]->values(),
                                cframes[attrib[i].cmp_frm]);
        }
        else { // PJSA 4-8-05: user input ctheta
          ele->setCompositeData2(1, 0, 0, coefData[attrib[i].cmp_attr]->values(),
                                 nodes, attrib[i].cmp_theta);
        }
      }
      else {
        LayInfo *li = layInfo[attrib[i].cmp_attr];
        if(li == 0) {
          fprintf(stderr," *** WARNING: Attribute found that refers to"
                         " nonexistant composite data: %d \n",
                         attrib[i].cmp_attr+1);
          continue;
        }
        // PJSA 3-31-05: set layer material properties if necessary
        for(int k=0; k<li->nLayers(); ++k) {
          int mid = li->getLayerMaterialId(k);
          // cerr << "k = " << k << ", mid = " << mid << endl;
          if(mid > -1) {
            LayMat *lmk = layMat[mid];
            li->setLayerMaterialProperties(k,lmk->data);
          }
        }
        // type is 3 for LAYC 2 for LAYN
        int type = 3 - li->getType();
        if(attrib[i].cmp_frm > -1) { // user input cframe
          ele->setCompositeData(type, li->nLayers(), li->values(), 0,
                                cframes[attrib[i].cmp_frm]);
        }
        else { // PJSA 4-8-05: user input ctheta
          ele->setCompositeData2(type, li->nLayers(), li->values(), 0,
                                 nodes, attrib[i].cmp_theta);
        }
      }
    }
  }
  if(structure_element_count > 0) {
    global_average_E /= double(structure_element_count);
    global_average_nu /= double(structure_element_count);
  }
  if(fluid_element_count > 0) {
    global_average_rhof /= double(fluid_element_count);
  }

  // setup beam element offsets
  for(vector<OffsetData>::iterator offIt = offsets.begin();
        	  offIt != offsets.end(); ++offIt) {
    for(int i = offIt->first; i <= offIt->last; ++i) {
      if(elemSet[i] == 0) {
	fprintf(stderr,
	   " *** WARNING: Setting up offset on non-existent element %d ***\n",
			i+1);
      } else
	elemSet[i]->setOffset(offIt->o);
    }
  }

  // Set new Material data in elements
  map<int,int>::iterator uIter = matUsage.begin();
  while(uIter != matUsage.end()) {
    int elemNum = uIter->first;
    int matNum = uIter->second;
    Element *ele = elemSet[ elemNum ];

    // Check if element exists
    if (ele == 0) {
      fprintf(stderr," *** WARNING: Material was found for"
                       " non existent element %d \n",attrib[i].nele+1);
      uIter++; // go to the next mapping
      continue;
    }
    map<int, Material *>::iterator matIter = materials.find(matNum);
    if(matIter == materials.end()) {
         fprintf(stderr," *** WARNING: Non Existent material (%d)"
                        " was assigned to element %d \n", matNum+1, elemNum+1);
    } else
      ele->setMaterial(matIter->second);
    uIter++;
  }
}

//----------------------------------------------------------------

int GeoSource::getNodes(CoordSet &nds)  {

  nds.nodeCopy(nodes.last(), nodes);

  // the routines calling this expect to know the total number of nodes
  // including internally created ones.
  return numNodes+numInternalNodes;
}

CoordSet& GeoSource::GetNodes() { return nodes; }

//----------------------------------------------------------------

/* If nElems and elemList exist, then we are getting Elems for the
   binary case.  In this case, we will not need to create a glToPck
   since this is only used when the decomposition from text input
   Also, the packed ElemSet will have sorted the phantom and real
   elements, so that the phantom elements will be at the end of the list */

int GeoSource::getElems(Elemset &packedEset, int nElems, int *elemList)
{
  SolverInfo sinfo = domain->solInfo();
  //ADDED FOR HEV PROBLEM, EC, 20070820
  if(sinfo.HEV) { packedEsetFluid = new Elemset(); nElemFluid = 0; }

  int iEle, numele;

  if(nElems) numele = nElems;
  else numele = elemSet.last();

  int packFlag = 0;
  if(!nElems) packFlag = 1;

  nElem = 0; // counting the phantom elements as well
  int nPhantoms = 0;

  // add real elements to list
  for(iEle = 0; iEle < numele; ++iEle) {
    Element *ele = elemSet[iEle];
    if(ele)  {
      //if(ele->isPhantomElement()==0) {//getProperty()  //original
      if(!ele->isPhantomElement() && (!ele->isHEVFluidElement() || (!sinfo.HEV))) {//getProperty()  //MODIFIED FOR HEV PROBLEM, EC, 20070820
        packedEset.elemadd(nElem, ele);
        packedEset[nElem]->setGlNum(iEle);
        if(packFlag)
          glToPckElems[iEle] = nElem;
        nElem++;
      }
      //ADDED FOR HEV PROBLEM, EC, 20070820
      //else if(!ele->isPhantomElement() && ele->isHEVFluidElement() && (sinfo.HEV)) {
      else if(!ele->isPhantomElement() && ele->isHEVFluidElement()) {
        packedEsetFluid->elemadd(nElemFluid, ele);
        (*packedEsetFluid)[nElemFluid]->setGlNum(iEle);
        nElemFluid++;
      }
      else
        nPhantoms++;
    }
  }

  if(!domain->getSowering())
    {
      //add the sommerElement from vector sommer -JF
      //cerr << "GeoSourceCore.C, adding sommer in packedEset" << endl;
      if (sinfo.ATDARBFlag!=-2.0) {
	for (int i = 0 ;i<domain->numSommer;i++) {
	  packedEset.elemadd(nElem,domain->sommer[i]);
	  packedEset[nElem]->setGlNum(numele+i);
	  if(packFlag)
	    glToPckElems[numele+i] = nElem;
	  nElem++;
	}
      }
    }

  int nRealElems = nElem;

  // set number of real elements
  packedEset.setEmax(nElem);
  packedEset.setNumPhantoms(nPhantoms);

  //ADDED FOR HEV PROBLEM, EC, 20070820
  if(sinfo.HEV) {
    packedEsetFluid->setEmax(nElemFluid);
    packedEsetFluid->setNumPhantoms(0);
  }

  // add phantom elements to list
  iEle = 0;
  int nPhants = 0;
  while (nPhants < nPhantoms)  {
    Element *ele = elemSet[iEle];
    if(ele)  {
      if(ele->isPhantomElement()) {//getProperty() == 0
        packedEset.elemadd(nElem, ele);
        packedEset[nElem]->setGlNum(iEle);
        if(packFlag)
          glToPckElems[iEle] = nElem;
        nElem++;
        nPhants++;
      }
    }
    iEle++;
  }
  return nRealElems;

}

/*
//identical to original getElems function; for generating top files
int GeoSource::getElemsTopHEV(Elemset &packedEset, int nElems, int *elemList)
{
  SolverInfo sinfo = domain->solInfo();
  //ADDED FOR HEV PROBLEM, EC, 20070820
  //if(sinfo.HEV) { packedEsetFluid = new Elemset(); nElemFluid = 0; }

  int iEle, numele;

  if(nElems) numele = nElems;
  else numele = elemSet.last();

  int packFlag = 0;
  if(!nElems) packFlag = 1;

  nElem = 0; // counting the phantom elements as well
  int nPhantoms = 0;

  // add real elements to list
  for(iEle = 0; iEle < numele; ++iEle) {
    Element *ele = elemSet[iEle];
    if(ele)  {
      if(ele->isPhantomElement()==0) {//getProperty()  //original
      //if(!ele->isPhantomElement() && !ele->isHEVFluidElement()) {//getProperty()  //MODIFIED FOR HEV PROBLEM, EC, 20070820
        packedEset.elemadd(nElem, ele);
        packedEset[nElem]->setGlNum(iEle);
        if(packFlag)
          glToPckElems[iEle] = nElem;
        nElem++;
      }
      //ADDED FOR HEV PROBLEM, EC, 20070820
      //else if(!ele->isPhantomElement() && ele->isHEVFluidElement()) {
        //packedEsetFluid->elemadd(nElemFluid, ele);
        //(*packedEsetFluid)[nElemFluid]->setGlNum(iEle);
        //nElemFluid++;
      //}
      else
        nPhantoms++;
    }
  }

  if(!domain->getSowering())
    {
      //add the sommerElement from vector sommer -JF
      //cerr << "GeoSourceCore.C, adding sommer in packedEset" << endl;
      if (sinfo.ATDARBFlag!=-2.0) {
	for (int i = 0 ;i<domain->numSommer;i++) {
	  packedEset.elemadd(nElem,domain->sommer[i]);
	  packedEset[nElem]->setGlNum(numele+i);
	  if(packFlag)
	    glToPckElems[numele+i] = nElem;
	  nElem++;
	}
      }
    }

  int nRealElems = nElem;

  // set number of real elements
  packedEset.setEmax(nElem);
  packedEset.setNumPhantoms(nPhantoms);

  //ADDED FOR HEV PROBLEM, EC, 20070820
  //if(sinfo.HEV) {
    //packedEsetFluid->setEmax(nElemFluid);
    //packedEsetFluid->setNumPhantoms(0);
  //}

  // add phantom elements to list
  iEle = 0;
  int nPhants = 0;
  while (nPhants < nPhantoms)  {
    Element *ele = elemSet[iEle];
    if(ele)  {
      if(ele->isPhantomElement()) {//getProperty() == 0
        packedEset.elemadd(nElem, ele);
        packedEset[nElem]->setGlNum(iEle);
        if(packFlag)
          glToPckElems[iEle] = nElem;
        nElem++;
        nPhants++;
      }
    }
    iEle++;
  }
  return nRealElems;

}
*/

int GeoSource::getNonMpcElems(Elemset &eset)
{
  int numele = elemSet.last();
  int mpccount = 0, elecount = 0;
  // add non-mpc elements to list
  for(int iEle = 0; iEle < numele; ++iEle) {
    Element *ele = elemSet[iEle];
    if(ele)  {
      if(!ele->isRigidMpcElement() && !ele->isMpcElement())  {
        eset.elemadd(elecount, ele);
        eset[nElem]->setGlNum(iEle);
        elecount++;
      }
      else mpccount++;
    }
    eset.setEmax(elecount);
  }
  return mpccount;
}

void GeoSource::setElemTypeMap()
{
  // build elem type to num nodes map
  elemTypeNumNodesMap = new int[89];
  for (int i = 0; i < 89; i++)
    elemTypeNumNodesMap[i] = -1;

  elemTypeNumNodesMap[1] = 2;
  elemTypeNumNodesMap[2] = 4;
  elemTypeNumNodesMap[3] = 4;
  elemTypeNumNodesMap[4] = 4;
  elemTypeNumNodesMap[6] = 2;
  elemTypeNumNodesMap[7] = 3;
  elemTypeNumNodesMap[8] = 3;
  elemTypeNumNodesMap[9] = 2;
  elemTypeNumNodesMap[10] = 4;
  elemTypeNumNodesMap[11] = 1;
  elemTypeNumNodesMap[17] = 8;
  elemTypeNumNodesMap[18] = 4;
  elemTypeNumNodesMap[19] = 3;
  elemTypeNumNodesMap[20] = 3;
  elemTypeNumNodesMap[21] = 2;
  elemTypeNumNodesMap[22] = 2;
  elemTypeNumNodesMap[23] = 4;
  elemTypeNumNodesMap[24] = 5;
  elemTypeNumNodesMap[25] = 10;
  elemTypeNumNodesMap[30] = 4;
  elemTypeNumNodesMap[31] = 4;
  elemTypeNumNodesMap[35] = 3;
  elemTypeNumNodesMap[36] = 3;
  elemTypeNumNodesMap[40] = 4;
  elemTypeNumNodesMap[41] = 4;
  elemTypeNumNodesMap[52] = 8;
  elemTypeNumNodesMap[60] = 4;
  elemTypeNumNodesMap[61] = 3;
  elemTypeNumNodesMap[88] = 4; //HB

  //return elemTypeNumNodesMap;
}

void GeoSource::setElementPressure(int elemNum, double pressure)
{
 prsflg = 1;

 if(elemSet[elemNum])
   elemSet[elemNum]->setPressure(pressure);
 else
   fprintf(stderr," *** WARNING: element %d does not exist \n", elemNum+1);
}

void GeoSource::setElementPreLoad(int elemNum, double preload)
{
 if(elemSet[elemNum])
   elemSet[elemNum]->setPreLoad(preload,prlflg);
 else
   fprintf(stderr," *** WARNING: element %d does not exist \n", elemNum+1);

}

void GeoSource::setConsistentPFlag()
{
 fprintf(stderr," ... Using Consistent Pressure Force ...\n");
 constpflg = 1;
}

void GeoSource::setConsistentQFlag()
{
 fprintf(stderr," ... Using Consistent Gravity Force ...\n");
 constqflg = 1;
}

int GeoSource::getTextDirichletBC(BCond *&bc)
{
  bc = textDBC;
  return numTextDirichlet;
}

int GeoSource::getTextNeumanBC(BCond *&bc)
{
  bc = textNBC;
  return numTextNeuman;
}

int GeoSource::getDirichletBC(BCond *&bc)
{
  bc = dbc;
  return numDirichlet;
}

int GeoSource::getDirichletBCFluid(BCond *&bc)
{
  if (dbcFluid) {
    bc = dbcFluid;
    return numDirichletFluid;
  }
  else {
   bc = 0;
   return 0;
  }
}

int GeoSource::getNeumanBC(BCond *&bc)
{
  bc = nbc;
  return numNeuman;
}

int GeoSource::getConvBC(BCond *&bc)
{
  bc = cvbc;
  return numConvBC;
}

int GeoSource::getRadBC(BCond *&bc)
{
  bc = rdbc;
  return numRadBC;
}

int GeoSource::getIDis(BCond *&bc)
{
  bc = iDis;
  return numIDis;
}

int GeoSource::getIDis6(BCond *&bc)
{
  bc = iDis6;
  return numIDis6;
}

int GeoSource::getIVel(BCond *&bc)
{
  bc = iVel;
  return numIVel;
}

int GeoSource::getITemp(BCond *&bc)
{
  bc = iTemp;
  return numITemp;
}

int GeoSource::getModalDamping(BCond *&damping)
{
  damping = modalDamping;
  return numDampedModes;
}

int GeoSource::getSurfaceDirichletBC(BCond *&bc)
{
  bc = surface_dbc;
  return numSurfaceDirichlet;
}

int GeoSource::getSurfaceNeumanBC(BCond *&bc)
{
  bc = surface_nbc;
  return numSurfaceNeuman;
}

int GeoSource::getSurfacePressure(BCond *&bc)
{
  bc = surface_pres;
  return numSurfacePressure;
}

void GeoSource::computeGlobalNumElements()
{
  // PJSA: note timing file will print nElem
  // need do a global max across all clusters
#ifdef DISTRIBUTED
  for(int i=0; i<numClusters; ++i)
    allNumClusElems[i] = structCom->globalMax(allNumClusElems[i]);  // find a better way!!
#endif
  nElemAllClusters = 0;
  for(int i=0; i<numClusters; ++i) nElemAllClusters += allNumClusElems[i];
  delete [] allNumClusElems; allNumClusElems = 0;
}

void GeoSource::setTextBC()
{
  // keep track of bc's defined in text file for addition to binary data
  numTextDirichlet = numDirichlet;
  if (numTextDirichlet)
    textDBC = dbc;
  numTextNeuman = numNeuman;
  if (numTextNeuman)
    textNBC = nbc;
}

void GeoSource::applyAuxData(int *cl2LocElem, int *cl2LocNode,
                             int minElemNum, int maxElemNum)
{
  // renumber elements in this subdomain
  for (int iElem = 0; iElem < nElem; iElem++)
    elemSet[iElem]->renum(cl2LocNode);

  // Set up element attributes

  // set material properties for element
  std::map<int, StructProp *> usedMat;
  for (int iAttr = minElemNum; iAttr < maxElemNum; iAttr++)  {

    int locElemNum = cl2LocElem[ attrib[iAttr].nele ];

    if (locElemNum >= 0)  {
      if (attrib[iAttr].attr >= 0)  {
        StructProp *prop;
        std::map<int, StructProp *>::iterator it = usedMat.find(attrib[iAttr].attr);
        if(it == usedMat.end())  {
          // prop = new StructProp(sProps[attrib[iAttr].attr]);
          prop = &sProps[attrib[iAttr].attr]; // PJSA: for helmsweep, helmCoef is updated in sProps[0] so don't copy!!
                                              // check with Thuan why he decided to make a copy of the material here
          usedMat[attrib[iAttr].attr] = prop;
        }
        else
          prop = it->second;

        elemSet[locElemNum]->setProp(prop);
      }
      else  {
        phantomFlag = 1;
        elemSet[locElemNum]->setProp(0);
      }
    }

    if(attrib[iAttr].cmp_attr >= 0) {
      if(coefData[attrib[iAttr].cmp_attr] != 0) {
        if(attrib[iAttr].cmp_frm > -1) { // user input cframe
          elemSet[locElemNum]->setCompositeData(1, 0, 0, coefData[attrib[iAttr].cmp_attr]->values(),
                                                cframes[attrib[iAttr].cmp_frm]);
        }
        else { // user input ctheta
          elemSet[locElemNum]->setCompositeData2(1, 0, 0, coefData[attrib[iAttr].cmp_attr]->values(),
                                                 nodes, attrib[iAttr].cmp_theta);
        }
      }
      else  {
        if (layInfo[attrib[iAttr].cmp_attr] == 0) {
          fprintf(stderr," *** WARNING: Attribute found that refers to"
                         " nonexistant composite data: %d"
                         " for elemNum %d \n",
                           attrib[iAttr].cmp_attr+1, iAttr);
          continue;
        }
      }

      // type is 3 for LAYC 2 for LAYN
      int type = 3-layInfo[attrib[iAttr].cmp_attr]->getType();
      if(attrib[iAttr].cmp_frm > -1) { // user input cframe
        elemSet[locElemNum]->setCompositeData(type, layInfo[attrib[iAttr].cmp_attr]->nLayers(), layInfo[attrib[iAttr].cmp_attr]->values(), 0,
                                              cframes[attrib[iAttr].cmp_frm]);
      }
      else { // user input ctheta
        elemSet[locElemNum]->setCompositeData2(type, layInfo[attrib[iAttr].cmp_attr]->nLayers(), layInfo[attrib[iAttr].cmp_attr]->values(), 0,
                                               nodes, attrib[iAttr].cmp_theta);
      }
    }
  }

  // set up element frames
  for(int iFrame = 0; iFrame < numEframes; iFrame++)  {

    // get local element number
    if(efd[iFrame].elnum >= maxElemNum)
      continue;

    int locElemNum =  cl2LocElem[ efd[iFrame].elnum ];

    // check if eframe is in this subdomain
    if(locElemNum >= 0)
      elemSet[locElemNum]->setFrame(&(efd[iFrame].frame));
  }
}

int GeoSource::getSubCtrl(BCond *glCtrlData, int numGlData,
               BCond *&locCtrlData, int glSub, int *&locToGlUserData)
{
  int i, iSub;

  // Form NodeToSub connectivity
  Connectivity *nodeToSub = subToNode->reverse();

  // count number of data in this subdomain
  int numLocCtrlData = 0;
  for (i = 0; i < numGlData; i++)  {
    for (iSub = 0; iSub < nodeToSub->num(glCtrlData[i].nnum); iSub++)
      if ( (*nodeToSub)[glCtrlData[i].nnum][iSub] == glSub )  {
        numLocCtrlData++;
        break;
      }
  }

  if (numLocCtrlData)  {

    // allocate arrays to number of data in this subdomain
    locCtrlData = new BCond[numLocCtrlData];
    locToGlUserData = new int[numLocCtrlData];

    // populate arrays
    int count = 0;
    for (i = 0; i < numGlData; i++)  {
      for (iSub = 0; iSub < nodeToSub->num(glCtrlData[i].nnum); iSub++)
        if ( (*nodeToSub)[glCtrlData[i].nnum][iSub] == glSub )  {

          locCtrlData[count] = glCtrlData[i];
          locToGlUserData[count] = i;

          // increment count
          count++;
        }
    }
  }

  return numLocCtrlData;
}

//---------------------------------------------------------------------

int GeoSource::getSubCtrl(BCond *glCtrlData, int numGlData,
                          BCond *&locCtrlData, int *gl2ClMap, int *cl2LocMap,
                          int *&locToGlUserData)
{
  int i;

  // count number of data in this subdomain
  int numLocCtrlData = 0;
  for (i = 0; i < numGlData; i++)  {
    if (glCtrlData[i].nnum <= maxGlobNode)  {
      int clusterNum = gl2ClMap[glCtrlData[i].nnum];
      if (clusterNum < maxClusNode && clusterNum >= 0)
        if (cl2LocMap[clusterNum] >= 0)
          numLocCtrlData++;
    }
  }

  if (numLocCtrlData)  {

    // allocate arrays to number of data in this subdomain
    locCtrlData = new BCond[numLocCtrlData];
    locToGlUserData = new int[numLocCtrlData];

    // populate arrays
    int count = 0;
    for (i = 0; i < numGlData; i++)  {

      if (glCtrlData[i].nnum <= maxGlobNode)  {
        int clusterNodeNum = gl2ClMap[glCtrlData[i].nnum];
        if (clusterNodeNum < maxClusNode && clusterNodeNum >= 0)  {
          if (cl2LocMap[clusterNodeNum] >= 0)  {

            locCtrlData[count] = glCtrlData[i];
            locToGlUserData[count] = i;

            // renumber to local node number
            locCtrlData[count].nnum = cl2LocMap[clusterNodeNum];

            // increment count
            count++;
          }
        }
      }
    }
  }
  return numLocCtrlData;
}

//--------------------------------------------------------------

void GeoSource::augmentBC(int numTextBC, BCond *textBC, BCond *&binBC,
                          int &numBinBC)  {

  BCond *tmpBC = 0;
  if (numBinBC)
    tmpBC = binBC;

  binBC = new BCond[numTextBC + numBinBC];

  int iBC;
  for (iBC = 0; iBC < numBinBC; iBC++)
    binBC[iBC] = tmpBC[iBC];

  for (iBC = 0; iBC < numTextBC; iBC++)
    binBC[numBinBC+iBC] = textBC[iBC];

  numBinBC += numTextBC;
}


//-------------------------------------------------------------------

int GeoSource::getBC(BCond *bc, int numBC, int *cl2LocNodeMap,
                     BCond *&subBC, int *gl2ClNodeMap){

  int numLocBC = 0;
  int bcNodeNum;

  // check if mapping exists for bc nodes
  int iBC;
  for (iBC = 0; iBC < numBC; iBC++)  {

    if (gl2ClNodeMap)  {
      if (bc[iBC].nnum <= maxGlobNode)
        bcNodeNum = gl2ClNodeMap[bc[iBC].nnum];
      else
        bcNodeNum = -1;
    }
    else
      bcNodeNum = bc[iBC].nnum;

    if (bcNodeNum >= 0  && bcNodeNum < maxClusNode)
      if (cl2LocNodeMap[ bcNodeNum ] >= 0)
        numLocBC++;
  }

  // set BC
  subBC = new BCond[numLocBC];

  int count = 0;
  for (iBC = 0; iBC < numBC; iBC++)  {

    if (gl2ClNodeMap)  {
      if (bc[iBC].nnum <= maxGlobNode)
        bcNodeNum = gl2ClNodeMap[bc[iBC].nnum];
      else
        bcNodeNum = -1;
    }
    else
      bcNodeNum = bc[iBC].nnum;

    if (bcNodeNum >= 0 && bcNodeNum < maxClusNode)  {

      if (cl2LocNodeMap[bcNodeNum] >= 0)  {

          subBC[count] = bc[iBC];

          // renumber to local node number
          subBC[count].nnum = cl2LocNodeMap[bcNodeNum];

          // increment count
          count++;
      }
    }
  }

  return numLocBC;
}

//--------------------------------------------------------------

void GeoSource::cleanAuxData()  {

  //layInfo.deleteArray();
  layInfo.restartArray();
  numLayInfo = 0;

  //coefData.deleteArray();
  coefData.restartArray();
  numCoefData = 0;

  numProps = 0;
  //sProps.deleteArray();
  //sProps.restartArray();

  na = 0;
  namax = 0;
  //attrib.deleteArray();
  //attrib.restartArray();

  numEframes = 0;

  cframes.restartArray();
  numCframes = 0;

  if (dbc)  {
    delete [] dbc;
    dbc = 0;
  }

  if (nbc)  {
    delete [] nbc;
    nbc = 0;
  }

  if (cvbc)  {
    delete [] cvbc;
    cvbc = 0;
  }

  if (iDis)  {
    delete [] iDis;
    iDis = 0;
  }

  if (iVel) {
    delete [] iVel;
    iVel = 0;
  }

  if (iTemp)  {
    delete [] iTemp;
    iTemp = 0;
  }
}

//--------------------------------------------------------------

int GeoSource::readRanges(BinFileHandler &file, int &numRanges,
                		int (*&ranges)[2])
{
  // read number of ranges
  file.read(&numRanges, 1);

  // read number of singles
  int numSingles;
  file.read(&numSingles, 1);

  // allocate memory for ranges
  ranges = new int[numRanges+numSingles][2];

  // read in ranges
  file.read(reinterpret_cast<int *>(ranges), 2*numRanges);

  // read the singles for this subdomain
  int *singles = new int[numSingles];
  file.read(singles, numSingles);

  // compute number of values in subdomain
  int numValues = 0;
  for (int iRange = 0; iRange < numRanges; iRange++)
    numValues += ranges[iRange][1] - ranges[iRange][0] + 1;

  numValues += numSingles;

  // convert singles to range format
  for (int iSingle = 0; iSingle < numSingles; iSingle++)  {

    ranges[numRanges+iSingle][0] = singles[iSingle];
    ranges[numRanges+iSingle][1] = singles[iSingle];
  }
  numRanges += numSingles;

  delete [] singles;
  return numValues;
}

//-----------------------------------------------------------------
void GeoSource::outputNodeVectors6(int fileNum, double (*xyz)[11],
	 			   int outputSize, double time)//DofSet::max_known_nonL_dof
{
  // 6 dof output should include node number (for IDISP6)
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time >= 0.0) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
    else
      filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
  }

  for(int inode = 0; inode < outputSize; inode++)  {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr,
             " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
             w, p, xyz[inode][0], w, p, xyz[inode][1],
             w, p, xyz[inode][2], w, p, xyz[inode][3],
             w, p, xyz[inode][4], w, p, xyz[inode][5]);
    else
      filePrint(oinfo[fileNum].filptr,
               " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
               w, p, xyz[inode][0], w, p, xyz[inode][1],
               w, p, xyz[inode][2], w, p, xyz[inode][3],
 	       w, p, xyz[inode][4], w, p, xyz[inode][5]);
  }
  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputNodeVectors6(int fileNum, DComplex (*xyz)[11],
                                   int outputSize, double time)//DofSet::max_known_nonL_dof
{
  // 6 dof output should include node number (for IDISP6)
  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
    default:
    case OutputInfo::realimag :
      // print real part or both real & imag parts for single node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E\n",w,p,time);
      }
      for(int inode = 0; inode < outputSize; inode++)  {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,
                  " % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E\n",
                  w, p, xyz[inode][0].real(), w, p, xyz[inode][1].real(),
                  w, p, xyz[inode][2].real(), w, p, xyz[inode][3].real(),
                  w, p, xyz[inode][4].real(), w, p, xyz[inode][5].real(),
                  w, p, xyz[inode][0].imag(), w, p, xyz[inode][1].imag(),
                  w, p, xyz[inode][2].imag(), w, p, xyz[inode][3].imag(),
                  w, p, xyz[inode][4].imag(), w, p, xyz[inode][5].imag());
        else
          filePrint(oinfo[fileNum].filptr,
                  " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                  w, p, xyz[inode][0].real(), w, p, xyz[inode][1].real(),
                  w, p, xyz[inode][2].real(), w, p, xyz[inode][3].real(),
                  w, p, xyz[inode][4].real(), w, p, xyz[inode][5].real());
      }
      // print imaginary part
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
        }
        for(int inode = 0; inode < outputSize; inode++)  {
          filePrint(oinfo[fileNum].filptr,
                    " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                    w, p, xyz[inode][0].imag(), w, p, xyz[inode][1].imag(),
                    w, p, xyz[inode][2].imag(), w, p, xyz[inode][3].imag(),
                    w, p, xyz[inode][4].imag(), w, p, xyz[inode][5].imag());
        }
      }
      break;
    case OutputInfo::modulusphase :
      // print modulus or modulus & phase for single node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E\n",w,p,time);
      }
      for(int inode = 0; inode < outputSize; inode++)  {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,
                  " % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E\n",
                  w, p, std::abs(xyz[inode][0]), w, p, std::abs(xyz[inode][1]),    //CRW
                  w, p, std::abs(xyz[inode][2]), w, p, std::abs(xyz[inode][3]),    //CRW
                  w, p, std::abs(xyz[inode][4]), w, p, std::abs(xyz[inode][5]),    //CRW
                  w, p, std::arg(xyz[inode][0]), w, p, arg(xyz[inode][1]),
                  w, p, std::arg(xyz[inode][2]), w, p, arg(xyz[inode][3]),
                  w, p, std::arg(xyz[inode][4]), w, p, arg(xyz[inode][5]));
        else
          filePrint(oinfo[fileNum].filptr,
                    " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                    w, p, std::abs(xyz[inode][0]), w, p, std::abs(xyz[inode][1]),    //CRW
                    w, p, std::abs(xyz[inode][2]), w, p, std::abs(xyz[inode][3]),    //CRW
                    w, p, std::abs(xyz[inode][4]), w, p, std::abs(xyz[inode][5]));    //CRW
      }
      // print phase
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E\n",w,p,time);
        }
        for(int inode = 0; inode < outputSize; inode++)  {
          filePrint(oinfo[fileNum].filptr,
                    " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                    w, p, arg(xyz[inode][0]), w, p, arg(xyz[inode][1]),
                    w, p, arg(xyz[inode][2]), w, p, arg(xyz[inode][3]),
                    w, p, arg(xyz[inode][4]), w, p, arg(xyz[inode][5]));
        }
      }
      break;
    case OutputInfo::animate :
      if(outputSize != 1) {
        double phi = 0;
        double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
        for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E\n",w,p,phi);
          for(int j = 0; j < outputSize; j++) {
            double proj[6];
            for(int k=0; k<6; ++k)
              proj[k] = std::abs(xyz[j][k])*cos(arg(xyz[j][k])-phi);    //CRW
            filePrint(oinfo[fileNum].filptr,
                      " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                      w, p, proj[0], w, p, proj[1], w, p, proj[2],
                      w, p, proj[3], w, p, proj[4], w, p, proj[5]);
          }
          phi += incr;
        }
      }
      else cerr << " *** WARNING: animate not supported for single-node output \n";
      break;
  }


  fflush(oinfo[fileNum].filptr);
}

//-----------------------------------------------------------------

// NOTE: This works only for 1 cluster

void GeoSource::outputNodeVectors(int fileNum, double (*glv)[11],
		                  int outputSize, double time)//DofSet::max_known_nonL_dof
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time >= 0.0) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
    else
      filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
  }

  int i;
  for(i = 0; i < outputSize; i++) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
              w,p,glv[i][0], w,p,glv[i][1], w,p,glv[i][2]);
    else
      filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                w,p,glv[i][0], w,p,glv[i][1], w,p,glv[i][2]);
  }
  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputNodeVectors(int fileNum, DComplex (*glv)[11],
                int outputSize, double time)//DofSet::max_known_nonL_dof
{
  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
    default:
    case OutputInfo::realimag :
      // print real part (or both real & imag in the case of 1 node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
      }
      for(i = 0; i < outputSize; i++)  {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr, " % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E \n",
                  w,p,glv[i][0].real(), w,p,glv[i][0].imag(), w,p,glv[i][1].real(),
                  w,p,glv[i][1].imag(), w,p,glv[i][2].real(), w,p,glv[i][2].imag());
        else
          filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                    w,p,glv[i][0].real(), w,p,glv[i][1].real(), w,p,glv[i][2].real());
      }
      // print imaginary part
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
        }
        for(i = 0; i < outputSize; i++)  {
          filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                    w,p,glv[i][0].imag(), w,p,glv[i][1].imag(), w,p,glv[i][2].imag());
        }
      }
      break;
    case OutputInfo::modulusphase :
      // print modulus or both modulus and phase in the case of 1 node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
      }
      for(i = 0; i < outputSize; i++)  {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr, " % *.*E % *.*E  % *.*E % *.*E  % *.*E % *.*E \n",
                  w,p,std::abs(glv[i][0]), w,p,std::abs(glv[i][1]), w,p,std::abs(glv[i][2]),    //CRW
                  w,p,arg(glv[i][0]), w,p,arg(glv[i][1]), w,p,arg(glv[i][2]));
        else
          filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                    w,p,std::abs(glv[i][0]), w,p,std::abs(glv[i][1]), w,p,std::abs(glv[i][2]));    //CRW
      }
      // print phase
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,time);
        }
        for(i = 0; i < outputSize; i++)  {
          filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                    w,p,arg(glv[i][0]), w,p,arg(glv[i][1]), w,p,arg(glv[i][2]));
        }
      }
      break;
    case OutputInfo::animate :
      if(outputSize != 1) {
        double phi = 0;
        double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
        for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,phi);
          for(int j = 0; j < outputSize; j++) {
            double proj[3];
            for(int k=0; k<3; ++k)
              proj[k] = std::abs(glv[j][k])*cos(arg(glv[j][k])-phi);    //CRW
            filePrint(oinfo[fileNum].filptr, " % *.*E % *.*E % *.*E\n",
                      w,p,proj[0], w,p,proj[1], w,p,proj[2]);
          }
          phi += incr;
        }
      }
      else cerr << " *** WARNING: animate not supported for single-node output \n";
      break;
  }

  fflush(oinfo[fileNum].filptr);
}

//------------------------------------------------------------

void GeoSource::outputNodeScalars(int fileNum, double *data,
				  int outputSize, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time >= 0.0) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
    else
      filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
  }

  for(int i = 0; i < outputSize; i++) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i]);
    else
      filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i]);
  }

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputNodeScalars(int fileNum, DComplex *data,
                                  int outputSize, double time)
{
  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
    default:
    case OutputInfo::realimag :
      // print real part
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E", w, p, data[i].real());
        else
          filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].real());
      }
      // print imaginary part
      if(time >= 0.0) {
        if(outputSize != 1) filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].imag());
        else
          filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].imag());
      }
      break;
    case OutputInfo::modulusphase :
      // print modulus
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E", w, p, std::abs(data[i]));     //CRW
        else
          filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, std::abs(data[i]));    //CRW
      }
      // print phase part
      if(time >= 0.0) {
        if(outputSize != 1) filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E\n", w, p, arg(data[i]));
        else
          filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, arg(data[i]));
      }
      break;
    case OutputInfo::animate :
      if(outputSize != 1) {
        double phi = 0;
        double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
        for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,phi);
          for(int j = 0; j < outputSize; j++) {
            double proj = std::abs(data[j])*cos(arg(data[j])-phi);    //CRW
            filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, proj);
          }
          phi += incr;
        }
      }
      else cerr << " *** WARNING: animate not supported for single-node output \n";
      break;
  }

  fflush(oinfo[fileNum].filptr);
}


//------------------------------------------------------------

void GeoSource::outputElemVectors(int fileNum, double *data,
                                  int outputSize, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time >= 0.0) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
    else
      filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
  }

  for(int i = 0; i < outputSize; i++) {
    if(outputSize == 1)
      fprintf(oinfo[fileNum].filptr," % *.*E % *.*E\n",
              w, p, data[2*i], w, p, data[2*i+1]);
    else
      filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                w, p, data[2*i], w, p, data[2*i+1]);
  }

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputElemVectors(int fileNum, DComplex *data,
                                  int outputSize, double time)
{
  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
    default:
    case OutputInfo::realimag :
      // output real part or both real & imag parts in the case of single node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E % *.*E  % *.*E % *.*E\n",
                  w, p, data[2*i].real(), w, p, data[2*i+1].real(),
                  w, p, data[2*i].imag(), w, p, data[2*i+1].imag());
        else
          filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                    w, p, data[2*i].real(), w, p, data[2*i+1].real());
      }
      // output imaginary part
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
        }
        for(i = 0; i < outputSize; i++)
          filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                    w, p, data[2*i].imag(), w, p, data[2*i+1].imag());
      }
      break;
   case OutputInfo::modulusphase :
      // output modulus part or both modulus & phase in the case of single node output
      if(time >= 0.0) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
        else
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
      }
      for(i = 0; i < outputSize; i++) {
        if(outputSize == 1)
          fprintf(oinfo[fileNum].filptr," % *.*E % *.*E  % *.*E % *.*E\n",
                  w, p, std::abs(data[2*i]), w, p, std::abs(data[2*i+1]),    //CRW
                  w, p, arg(data[2*i]), w, p, arg(data[2*i+1]));
        else
          filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                    w, p, std::abs(data[2*i]), w, p, std::abs(data[2*i+1]));    //CRW
      }
      // output phase part
      if(outputSize != 1) {
        if(time >= 0.0) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
        }
        for(i = 0; i < outputSize; i++)
          filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                    w, p, arg(data[2*i]), w, p, arg(data[2*i+1]));
      }
      break;
    case OutputInfo::animate :
      if(outputSize != 1) {
        double phi = 0;
        double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
        for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
          filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,phi);
          for(int j = 0; j < outputSize; j++) {
            double proj[2];
            for(int k=0; k<2; ++k)
              proj[k] = std::abs(data[2*j+k])*cos(arg(data[2*j+k])-phi);    //CRW
            filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
                      w, p, proj[0], w, p, proj[1]);
          }
          phi += incr;
        }
      }
      else cerr << " *** WARNING: animate not supported for single-node output \n";
      break;
  }

  fflush(oinfo[fileNum].filptr);
}

//-----------------------------------------------------------------

// now same function exists not reading from file
void GeoSource::getTextDecomp(bool sowering)
{
  MatrixTimers &mt = domain->getTimers();
  startTimerMemory(mt.readDecomp, mt.memoryDecomp);

  // Get Decomposition File pointer
  FILE *f = cinfo->decPtr;

  if(f == 0)
    f = fopen("DECOMPOSITION","r");

    if(f == 0) {
      fprintf(stderr," **************************************************\n");
      fprintf(stderr," *** ERROR: DECOMPOSITION file does not exist   ***\n");
      fprintf(stderr," ***        Please provide a DECOMPOSITION file ***\n");
      fprintf(stderr," **************************************************\n");
      exit(-1);
    }

  // get decomposition
  if(subToElem) {
    fprintf(stderr," ... Already read Decomposition\n");
    return;
  }
  mt.memorySubToElem -= memoryUsed();
  // Allocate memory for subdomain to element connectivity
  int *connect = new int[nElem];

  // Get the number of subdomains in Decomposition file
  int numSub;
  int error = fscanf(f,"%d",&numSub);

  // Decomposition file error checking
  if(error == 0) {
    char s1[14],s2[40],s3[4],s4[40];
    int error = fscanf(f,"%s%s%s%s",s1,s2,s3,s4);

    // Get the number of subdomains in Decomposition file
    error = fscanf(f,"%d",&numSub);
  }

  int *cx = new int[numSub+1];

  int curEle = 0;
  int isub;
  for (isub = 0; isub < numSub; ++isub) {
    int nele;
    fscanf(f,"%d",&nele);
    cx[isub] = curEle;
    if(curEle + nele > numElem()) {
      fprintf(stderr," *** ERROR: This decomposition contains more elements "
                     "than the Original mesh:\n");
      fprintf(stderr," *** %d vs %d\n", curEle + nele, nElem);
      exit(1);
    }

    int iele;
    for (iele = 0; iele < nele; ++iele) {
      fscanf(f,"%d",connect+curEle);
      connect[curEle] -= 01;
      curEle++;
    }
  }

  cx[numSub] = curEle;
  // PHIL FIX THIS
  /*int iele = 0;
  int cEle = 0;
  for(isub = 0; isub < numSub; ++isub) {
    for(; iele < cx[isub+1]; iele++)
      if(glToPckElems.find(connect[iele]) != glToPckElems.end())
        connect[cEle++] =  connect[iele];
    cx[isub+1] = cEle;
    }*/

  subToElem = new Connectivity(numSub,cx,connect);

  subToElem->renumberTargets(glToPckElems);

  mt.memorySubToElem += memoryUsed();

  stopTimerMemory(mt.readDecomp, mt.memoryDecomp);
#ifdef DISTRIBUTED
  if(!binaryInput) {
    int* ptr = new int[numSub+1];
    int* target = new int[numSub];
    ptr[0] = 0;
    for(int i=0;i<numSub; i++){
       ptr[i+1] = ptr[i]+1;
       target[i] = 0;
    }
    subToClus = new Connectivity(numSub,ptr,target);
    numClusters = 1;
    clusToSub = subToClus->reverse();
    numClusNodes = nGlobNodes;
    numClusElems = nElem; //HB: not sure this is be always correct (i.e. phantoms els ...)
  }
  //setNumNodalOutput();// -JFD
#endif
}

void GeoSource::setNumNodalOutput()
{
  // PJSAX: for single node output
  for(int iInfo = 0; iInfo < numOutInfo; iInfo++)
    if(oinfo[iInfo].nodeNumber != -1)
      numNodalOutput++;
  if(numNodalOutput)  {
    outputNodes = new int[numNodalOutput];
    outNodeIndex = new int[numNodalOutput];
    int count = 0;
    for(int iInfo = 0; iInfo < numOutInfo; iInfo++)
      if(oinfo[iInfo].nodeNumber != -1)  {
        outputNodes[count] = oinfo[iInfo].nodeNumber;
        outNodeIndex[count] = iInfo;
        count++;
      }
  }
}

//----------------------------------------------------------------------

int GeoSource::setDirichlet(int _numDirichlet, BCond *_dbc)
{
  if(domain->solInfo().fetiInfo.dmpc) {
    domain->addDirichletLMPCs(_numDirichlet, _dbc);
    return 0;
  }
  /*
  // MODIFIED FOR HEV PROBLEM, EC, 20070820
  if(domain->solInfo().HEV) {
      int i;
      int idbcFluid = 0;
      int idbcOthers = 0;

      BCond *nd_tmp = new BCond[_numDirichlet];
      BCond *ndFluid_tmp = new BCond[_numDirichlet];

      for(i=0;i<_numDirichlet;++i) {
        fprintf(stderr," ... dofnum of BC %d is %d ...\n",i+1,_dbc[i].dofnum);
        if(_dbc[i].dofnum == 10) {
          fprintf(stderr," ... Found Fluid Dirichlet BC ...\n");
          ndFluid_tmp[idbcFluid] = _dbc[i];
          idbcFluid++;
        }
        else {
          nd_tmp[idbcOthers] = _dbc[i];
          idbcOthers++;
        }
      }
      // Allocate memory for correct number of dbc
      BCond *nd = new BCond[numDirichlet+idbcOthers];
      BCond *ndFluid = new BCond[numDirichletFluid+idbcFluid];


      // copy old dbcFluid
      if (dbcFluid) {
        for(i = 0; i<numDirichletFluid; ++i)
          ndFluid[i] = dbcFluid[i];

        // copy new dbcFluid
        for(i = 0; i<idbcFluid; ++i)
          ndFluid[i+numDirichletFluid] = ndFluid_tmp[i];

        // set correct number of dbcFluid
        numDirichletFluid += idbcFluid;

        // delete old array of dbcFluid
        delete [] dbcFluid;
        delete [] ndFluid_tmp;

        // set new pointer to correct number of dbcFluid
        dbcFluid = ndFluid;
      }
      else {
        for(i=0;i<idbcFluid;++i) {
          ndFluid[i] = ndFluid_tmp[i];
        }
        delete [] ndFluid_tmp;
        numDirichlet = idbcFluid;
        dbcFluid = ndFluid;
      }

      if (dbc) {
        for(i = 0; i < numDirichlet; ++i)
           nd[i] = dbc[i];

        // copy new dbc
        for(i = 0; i<idbcOthers; ++i)
          nd[i+numDirichlet] = nd_tmp[i];

        // set correct number of dbc
        numDirichlet += idbcOthers;

        // delete old array of dbc
        delete [] dbc;
        delete [] nd_tmp;

        // set new pointer to correct number of dbc
        dbc = nd;
      }
      else {
        for(i=0;i<idbcOthers;++i) {
          nd[i] = nd_tmp[i];
        }
        delete [] nd_tmp;
        numDirichlet = idbcOthers;
        dbc          = nd;
      }
  }
  else {
  */

  if(dbc) {

    // Allocate memory for correct number of dbc
    BCond *nd = new BCond[numDirichlet+_numDirichlet];

    // copy old dbc
    int i;
    for(i = 0; i < numDirichlet; ++i)
       nd[i] = dbc[i];

    // copy new dbc
    for(i = 0; i<_numDirichlet; ++i)
      nd[i+numDirichlet] = _dbc[i];

    // set correct number of dbc
    numDirichlet += _numDirichlet;

    // delete old array of dbc
    delete [] dbc;

    // set new pointer to correct number of dbc
    dbc = nd;

  }

  else {
    numDirichlet = _numDirichlet;
    dbc          = _dbc;
  }

  return 0;
}


//-------------------------------------------------------------------
void GeoSource::convertHEVDirToHelmDir()
{
  if(dbcFluid) {

    // Allocate memory for correct number of dbc
    BCond *nd = new BCond[numDirichlet+numDirichletFluid];

    // copy old dbc
    int i;
    for(i = 0; i < numDirichlet; ++i)
       nd[i] = dbc[i];

    // copy new dbc and modify the dofs
    for(i = 0; i<numDirichletFluid; ++i) {
      nd[i+numDirichlet] = dbcFluid[i];
      nd[i+numDirichlet].dofnum = 8-1;
    }

    // set correct number of dbc
    numDirichlet += numDirichletFluid;

    // delete old array of dbc
    if (dbc) delete [] dbc;
    delete[] dbcFluid;

    // set new pointer to correct number of dbc
    dbc = nd;
    numDirichletFluid = 0;
  }
}

//-------------------------------------------------------------------
int GeoSource::setDirichletFluid(int _numDirichletFluid, BCond *_dbcFluid)
{
  if(dbcFluid) {

    // Allocate memory for correct number of dbcFluid
    BCond *ndFluid = new BCond[numDirichletFluid+_numDirichletFluid];

    // copy old dbcFluid
    int i;
    for(i = 0; i < numDirichletFluid; ++i)
       ndFluid[i] = dbcFluid[i];

    // copy new dbcFluid
    for(i = 0; i<_numDirichletFluid; ++i)
      ndFluid[i+numDirichletFluid] = _dbcFluid[i];

    // set correct number of dbcFluid
    numDirichletFluid += _numDirichletFluid;

    // delete old array of dbcFluid
    delete [] dbcFluid;

    // set new pointer to correct number of dbcFluid
    dbcFluid = ndFluid;

  }

  else {
    numDirichletFluid = _numDirichletFluid;
    dbcFluid          = _dbcFluid;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setNeuman(int _numNeuman, BCond *_nbc)
{
  if (nbc) {

    // Allocate memory for correct number of nbc
    BCond *nd = new BCond[numNeuman+_numNeuman];

    // copy old nbc
    int i;
    for (i = 0; i < numNeuman; ++i)
      nd[i] = nbc[i];

    // copy new nbc
    for (i = 0; i < _numNeuman; ++i)
      nd[i+numNeuman] = _nbc[i];

    // set correct number of nbc
    numNeuman += _numNeuman;

    // delete old array of nbc
    delete [] nbc;

    // set new pointer to correct number of nbc
    nbc = nd;

  }
  else  {
    numNeuman = _numNeuman;
    nbc       = _nbc;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setConvBC(int _numConvBC, BCond *_cvbc)
{
  if (cvbc) {

    // Allocate memory for correct number of cvbc
    BCond *nd = new BCond[numConvBC+_numConvBC];

    // copy old cvbc
    int i;
    for (i = 0; i < numConvBC; ++i)
      nd[i] = cvbc[i];

    // copy new cvbc
    for (i = 0; i < _numConvBC; ++i)
      nd[i+numConvBC] = _cvbc[i];

    // set correct number of cvbc
    numConvBC += _numConvBC;

    // delete old array of cvbc
    delete [] cvbc;

    // set new pointer to correct number of cvbc
    cvbc = nd;

  }
  else  {
    numConvBC = _numConvBC;
    cvbc       = _cvbc;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setRadBC(int _numRadBC, BCond *_rdbc)
{
  if (rdbc) {

    // Allocate memory for correct number of rdbc
    BCond *nd = new BCond[numRadBC+_numRadBC];

    // copy old rdbc
    int i;
    for (i = 0; i < numRadBC; ++i)
      nd[i] = rdbc[i];

    // copy new rdbc
    for (i = 0; i < _numRadBC; ++i)
      nd[i+numRadBC] = _rdbc[i];

    // set correct number of rdbc
    numRadBC += _numRadBC;

    // delete old array of rdbc
    delete [] rdbc;

    // set new pointer to correct number of rdbc
    rdbc = nd;

  }
  else  {
    numRadBC = _numRadBC;
    rdbc       = _rdbc;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIDis6(int _numIDis6, BCond *_iDis6)
{
/*
  if (iDis6) {

    // Allocate memory for correct number of iDis6
    BCond *nd = new BCond[numIDis6+_numIDis6];

    // copy old iDis6
    int i;
    for (i = 0; i < numIDis6; ++i)
      nd[i] = iDis6[i];

    // copy new iDis6
    for ( i = 0; i < _numIDis6; ++i)
      nd[i+numIDis6] = _iDis6[i];

    // set correct number of iDis6
    numIDis6 += _numIDis6;

    // delete old array of iDis6
    delete [] iDis6;

    // set new pointer to correct number of iDis6
    iDis6 = nd;

  }
  else  {
*/

    numIDis6 = _numIDis6;
    iDis6    = _iDis6;
//  }
  return 0;
}

//-------------------------------------------------------------------
int GeoSource::setPitaIDis6(int n, BCond *i, int numTSPitaIDis6_)
{
  // numTSPitaIDis6 corresponds to the number of time-slices and also
	// to the number of displacement vectors read from the input file.
  // numPitaIDis6 corresponds to the number of boundary conditions for
	// each vector, that is 6 * #dof.
  // These values will be used to check that the all necessary seed
	// values have been specified
  numPitaIDis6   = n / numTSPitaIDis6_;
  numTSPitaIDis6 = numTSPitaIDis6_;
  PitaIDis6 = i;
  newStep0 = true;
  return 0;
}

//-------------------------------------------------------------------
int GeoSource::setPitaIVel6(int n, BCond *i, int numTSPitaIVel6_)
{
  // Same principle as setPitaIDis6()
  numPitaIVel6   = n / numTSPitaIVel6_;
  numTSPitaIVel6 = numTSPitaIVel6_;
  PitaIVel6 = i;
  newStep0 = true;
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIDis(int _numIDis, BCond *_iDis)
{
  if (iDis) {
    // Allocate memory for correct number of iDis
    BCond *nd = new BCond[numIDis+_numIDis];

    // copy old iDis
    for (int i = 0; i < numIDis; ++i)
      nd[i] = iDis[i];

    // copy new iDis
    for (int i = 0; i < _numIDis; ++i)
      nd[i+numIDis] = _iDis[i];

    // set correct number of iDis
    numIDis += _numIDis;

    // delete old array of iDis
    delete [] iDis;

    // set new pointer to correct number of iDis
    iDis = nd;
  }
  else {
    numIDis = _numIDis;
    iDis    = _iDis;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIVel(int _numIVel, BCond *_iVel)
{
  if (iVel) {

    // Allocate memory for correct number of iVel
    BCond *nd = new BCond[numIVel+_numIVel];

    // copy old iVel
    int i;
    for (i = 0; i < numIVel; ++i)
      nd[i] = iVel[i];

    // copy new iVel
    for ( i = 0; i < _numIVel; ++i)
      nd[i+numIVel] = _iVel[i];

    // set correct number of iVel
    numIVel += _numIVel;

    // delete old array of iVel
    delete [] iVel;

    // set new pointer to correct number of iVel
    iVel = nd;

  }
  else  {
    numIVel = _numIVel;
    iVel    = _iVel;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setITemp(int _numITemp, BCond *_iTemp)
{
  if (iTemp) {

    // Allocate memory for correct number of iTemp
    BCond *nd = new BCond[numITemp+_numITemp];

    // copy old iTemp
    int i;
    for (i = 0; i < numITemp; ++i)
      nd[i] = iTemp[i];

    // copy new iTemp
    for ( i = 0; i < _numITemp; ++i)
      nd[i+numITemp] = _iTemp[i];

    // set correct number of iTemp
    numITemp += _numITemp;

    // delete old array of iTemp
    delete [] iTemp;

    // set new pointer to correct number of iTemp
    iTemp = nd;

  }
  else  {
    numITemp = _numITemp;
    iTemp = _iTemp;
  }
  return 0;
}

//----------------------------------------------------------------

int GeoSource::setModalDamping(int _numDampedModes, BCond *_modalDamping)
{
  if(modalDamping) {

    // Allocate memory for correct number of damped modes
    BCond *nd = new BCond[numDampedModes + _numDampedModes];

    // copy old modalDamping
    int i;
    for(i = 0; i < numDampedModes; ++i)
      nd[i] = modalDamping[i];

    for(i = 0; i < _numDampedModes; ++i)
      nd[i+numDampedModes] = _modalDamping[i];

    numDampedModes += _numDampedModes;

    delete[] modalDamping;
    modalDamping = nd;
  }
  else {

    numDampedModes = _numDampedModes;
    modalDamping = _modalDamping;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfaceDirichlet(int _numSurfaceDirichlet, BCond *_surface_dbc)
{
  if(surface_dbc) {

    // Allocate memory for correct number of dbc
    BCond *nd = new BCond[numSurfaceDirichlet+_numSurfaceDirichlet];

    // copy old dbc
    int i;
    for(i = 0; i < numSurfaceDirichlet; ++i)
       nd[i] = surface_dbc[i];

    // copy new dbc
    for(i = 0; i<_numSurfaceDirichlet; ++i)
      nd[i+numSurfaceDirichlet] = _surface_dbc[i];

    // set correct number of dbc
    numSurfaceDirichlet += _numSurfaceDirichlet;

    // delete old array of dbc
    delete [] surface_dbc;

    // set new pointer to correct number of dbc
    surface_dbc = nd;

  }

  else {
    numSurfaceDirichlet = _numSurfaceDirichlet;
    surface_dbc          = _surface_dbc;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfaceNeuman(int _numSurfaceNeuman, BCond *_surface_nbc)
{
  if(surface_nbc) {

    // Allocate memory for correct number of nbc
    BCond *nd = new BCond[numSurfaceNeuman+_numSurfaceNeuman];

    // copy old nbc
    int i;
    for(i = 0; i < numSurfaceNeuman; ++i)
       nd[i] = surface_nbc[i];

    // copy new nbc
    for(i = 0; i<_numSurfaceNeuman; ++i)
      nd[i+numSurfaceNeuman] = _surface_nbc[i];

    // set correct number of nbc
    numSurfaceNeuman += _numSurfaceNeuman;

    // delete old array of nbc
    delete [] surface_nbc;

    // set new pointer to correct number of nbc
    surface_nbc = nd;

  }

  else {
    numSurfaceNeuman = _numSurfaceNeuman;
    surface_nbc          = _surface_nbc;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfacePressure(int _numSurfacePressure, BCond *_surface_pres)
{
  if(surface_pres) {

    // Allocate memory for correct number of pres
    BCond *nd = new BCond[numSurfacePressure+_numSurfacePressure];

    // copy old pres
    int i;
    for(i = 0; i < numSurfacePressure; ++i)
       nd[i] = surface_pres[i];

    // copy new pres
    for(i = 0; i<_numSurfacePressure; ++i)
      nd[i+numSurfacePressure] = _surface_pres[i];

    // set correct number of pres
    numSurfacePressure += _numSurfacePressure;

    // delete old array of pres
    delete [] surface_pres;

    // set new pointer to correct number of pres
    surface_pres = nd;

  }

  else {
    numSurfacePressure = _numSurfacePressure;
    surface_pres       = _surface_pres;
  }

  return 0;
}

//-------------------------------------------------------------------


void GeoSource::readMatchInfo(BinFileHandler &matchFile,
	int (*matchRanges)[2], int numMatchRanges, int subNum,
	int *clusToLocElem)  {

  // get TOC
  BinFileHandler::OffType tocLoc;
  matchFile.read(&tocLoc, 1);
  matchFile.seek(tocLoc);
  int numClusMatches;
  matchFile.read(&numClusMatches, 1);
  BinFileHandler::OffType dataStart;
  matchFile.read(&dataStart, 1);

  // get matches according to ranges
  for (int iR = 0; iR < numMatchRanges; iR++)  {

    // get number of data in range
    int nMatch = matchRanges[iR][1] - matchRanges[iR][0] + 1;

    // seek to correct file position
    matchFile.seek(dataStart + matchRanges[iR][0] * ( 2*sizeof(int) + 5*sizeof(double) ));

    for (int iMatch = 0; iMatch < nMatch; iMatch++)  {

      int tmpData[2];  // cluster elem #, global match file position
      matchFile.read(tmpData, 2);

      matchData[subNum][iMatch].elemNum = clusToLocElem[tmpData[0]];
      matchData[subNum][iMatch].glPos = tmpData[1];

      // read natual coordinates of match
      double coords[2];
      matchFile.read(coords, 2);

      matchData[subNum][iMatch].xi = coords[0];
      matchData[subNum][iMatch].eta = coords[1];

      // read gap vector
      matchFile.read(gapVec[subNum][iMatch], 3);
    }
  }
}

//-------------------------------------------------------------------

int GeoSource::getCPUMap(FILE *f, int numSub)
{
  int totSub = numSub;
  int numCPU;
  int *connect = new int[totSub];

#ifdef DISTRIBUTED
  if(f == 0) { // Trivial map

     numCPU = structCom->numCPUs();
     filePrint(stderr, " ... Making Trivial CPU Map, numCPU = %d ... \n", numCPU);
     int *cx  = new int[numCPU+1];
     subToCPU = new int[totSub];

     int curSub = 0;
     int icpu;
     int nSubPerCPU = totSub/numCPU;
     int remain = totSub%numCPU;
     for (icpu = 0; icpu < numCPU; ++icpu) {
       cx[icpu] = curSub;
       int iSub;
       int subForCPU = (icpu < remain) ? nSubPerCPU+1 : nSubPerCPU;
       for(iSub = 0; iSub < subForCPU; ++iSub) {
         connect[curSub] = curSub;
         subToCPU[curSub] = icpu;
         curSub++;
       }
     }
     cx[numCPU] = curSub;
     cpuToSub = new Connectivity(numCPU,cx,connect);
     // cpuToSub->print();

     Connectivity *subDomainToCPU = cpuToSub->reverse();
     cpuToCPU = cpuToSub->transcon(subDomainToCPU);
     delete subDomainToCPU;
  }
  else {
#endif
  int error = fscanf(f,"%d",&numCPU);
  filePrint(stderr, " ... Reading CPU Map from file %s, numCPU = %d ... \n", mapName, numCPU);
  if (error == 0) {
#ifdef DISTRIBUTED
    filePrint(stderr," *** ERROR: First line of CPU Map file is incorrect.\n");
    exit(-1);
#endif
  }

#ifdef DISTRIBUTED
  if(numCPU != structCom->numCPUs()) {
    fprintf(stderr, " *** ERROR: CPUMAP file is for %d MPI processes\n", numCPU);
    exit(-1);
  }
#endif

  int *cx  = new int[numCPU+1];
  subToCPU = new int[totSub];

  int curSub = 0;
  int icpu;
  for (icpu = 0; icpu < numCPU; ++icpu) {
    int nsub;
    fscanf(f,"%d",&nsub);
    cx[icpu] = curSub;
    if (curSub + nsub > totSub) {
      fprintf(stderr, "Exiting\n");
#ifdef USE_MPI
      if (structCom->myID() == 0)
#endif
        fprintf(stderr," *** ERROR: This decomposition contains more subdomains"
                       " than the Original mesh %d vs %d \n", curSub + nsub, totSub);
      fflush(stderr);
      exit(1);
    }
    int isub;
    for (isub=0; isub < nsub; ++isub) {
      fscanf(f,"%d",connect+curSub);
      connect[curSub] -= 1;
      subToCPU[connect[curSub]] = icpu;
      curSub++;
    }
  }
  cx[numCPU] = curSub;
  cpuToSub = new Connectivity(numCPU,cx,connect);

  Connectivity *subDomainToCPU = cpuToSub->reverse();
  cpuToCPU = cpuToSub->transcon(subDomainToCPU);
  delete subDomainToCPU;

#ifdef DISTRIBUTED
  }
#endif

  if(domain->solInfo().aeroFlag >= 0) {
    int numLocSub = 0;
#ifdef USE_MPI
    int myID = structCom->myID();
#else
    int myID = 0;
#endif
    numLocSub = cpuToSub->num(myID);

    // allocate for match data arrays
    typedef double (*gVec)[3];
    numMatchData = new int[numLocSub];
    matchData = new MatchData *[numLocSub];
    numGapVecs = new int[numLocSub];
    gapVec = new gVec[numLocSub];
  }

  return numCPU;
}

//-----------------------------------------------------------------------

void GeoSource::setMatchArrays(int numLocSub)  {

  // allocate for match data arrays
  typedef double (*gVec)[3];
  numMatchData = new int[numLocSub];
  matchData = new MatchData *[numLocSub];
  numGapVecs = new int[numLocSub];
  gapVec = new gVec[numLocSub];
}

//-----------------------------------------------------------------------
void GeoSource::addOutput(OutputInfo &outputInfo)  {

  oinfo[numOutInfo++] = outputInfo;
  if(outputInfo.type == OutputInfo::Farfield) domain->solInfo().farfield = true;
}

//-----------------------------------------------------------------------
std::pair<int, int>
GeoSource::getTimeSliceOutputFileIndices(int timeSliceRank) {
  std::map<int, std::pair<int, int> >::const_iterator it = timeSliceOutputFiles.find(timeSliceRank);
  return it != timeSliceOutputFiles.end() ? it->second : std::pair<int, int>(0, 0);
}

//--------------------------------------------------------------------
void GeoSource::openOutputFilesForPita(int sliceRank)
{
  std::pair<int, int> indices = getTimeSliceOutputFileIndices(sliceRank);
  for(int iInfo = indices.first; iInfo < indices.second; ++iInfo) {
   if (oinfo[iInfo].interval != 0) {
      char *fileName = oinfo[iInfo].filename;
      if (strlen(cinfo->outputExt) != 0) {
        int len1 = strlen(fileName);
        int len2 = strlen(cinfo->outputExt);
        char *nfn = new char[len1+len2+1];
        strcpy(nfn, fileName);
        strcat(nfn,cinfo->outputExt);
        fileName = nfn;
      }

      if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 )
        fprintf(stderr," *** ERROR: Cannot open %s\n",oinfo[iInfo].filename );
      outputHeader(iInfo);
      fflush(oinfo[iInfo].filptr);
    }
  }
}


//--------------------------------------------------------------------
void GeoSource::openOutputFiles(int *outNodes, int *outIndex, int numOuts)
{
  int iInfo;

  if(numOuts == 0) { // open all output files and write their corresponding TOPDOM/DEC header
    for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
      if(oinfo[iInfo].interval != 0) {
        char *fileName = oinfo[iInfo].filename;
        if (strlen(cinfo->outputExt) != 0) {
          int len1 = strlen(fileName);
          int len2 = strlen(cinfo->outputExt);
          char *nfn = new char[len1+len2+1];
          strcpy(nfn, fileName);
          strcat(nfn,cinfo->outputExt);
          fileName = nfn;
        }

        if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 )
          fprintf(stderr," *** ERROR: Cannot open %s\n",oinfo[iInfo].filename );
        outputHeader(iInfo);
        fflush(oinfo[iInfo].filptr);
      }
    }
  }
  else { // open selected output files
    for(int iOut = 0; iOut < numOuts; iOut++)  {
      iInfo = outIndex[iOut];
      if(oinfo[iInfo].interval != 0) {
        char *fileName = oinfo[iInfo].filename;
        if(strlen(cinfo->outputExt) != 0) {
          int len1 = strlen(fileName);
          int len2 = strlen(cinfo->outputExt);
          char *nfn = new char[len1+len2+1];
          strcpy(nfn, fileName);
          strcat(nfn,cinfo->outputExt);
          fileName = nfn;
        }

        if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 )
          fprintf(stderr," *** ERROR: Cannot open %s\n",oinfo[iInfo].filename );
        outputHeader(iInfo);
        fflush(oinfo[iInfo].filptr);
      }
    }
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFiles()
{
  for(int i = 0; i < numOutInfo; ++i) {
    this->closeOutputFileImpl(i);
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFilesForPita(int sliceRank)
{
  std::pair<int, int> indices = getTimeSliceOutputFileIndices(sliceRank);
  for (int i = indices.first; i < indices.second; ++i) {
    this->closeOutputFileImpl(i);
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFileImpl(int fileIndex)
{
  if((oinfo[fileIndex].interval != 0) && oinfo[fileIndex].filptr) {
    fclose(oinfo[fileIndex].filptr);
    oinfo[fileIndex].filptr = 0;
  }
}

//--------------------------------------------------------------------

void GeoSource::outputHeader(int fileNumber)
{
  // only one node is requested for output,
  if(oinfo[fileNumber].nodeNumber != -1) {
    fprintf(oinfo[fileNumber].filptr, "# node %d\n", oinfo[fileNumber].nodeNumber+1);
    return;
  }

  // get data description
  int headerSize = 200;
  char *headDescrip = new char[headerSize];   // description of data
  getHeaderDescription(headDescrip, fileNumber);
  fprintf(oinfo[fileNumber].filptr, headDescrip);

  delete [] headDescrip;
}

//---------------------------------------------------------------

void GeoSource::setControl(char *_checkfile, char *_nodeSetName,
                           char *_elemSetName, char *_bcondSetName)
{
/*
  cinfo->checkfile   = _checkfile;
  cinfo->nodeSetName = _nodeSetName;
  cinfo->elemSetName = _elemSetName;
  cinfo->bcondSetName = _bcondSetName;
*/
  // PJSA: limit these 4 strings to 31 characters + terminating null-character for xpost compatability
  cinfo->checkfile = new char [32];
  for(int i=0; i<32; ++i) {
    if(i==31) cinfo->checkfile[i] = '\0';
    else cinfo->checkfile[i] = _checkfile[i];
    if(_checkfile[i] == '\0') break;
  }
  cinfo->nodeSetName = new char [32];
  for(int i=0; i<32; ++i) {
    if(i==31) cinfo->nodeSetName[i] = '\0';
    else cinfo->nodeSetName[i] = _nodeSetName[i];
    if(_nodeSetName[i] == '\0') break;
  }
  cinfo->elemSetName = new char [32];
  for(int i=0; i<32; ++i) {
    if(i==31) cinfo->elemSetName[i] = '\0';
    else cinfo->elemSetName[i] = _elemSetName[i];
     if(_elemSetName[i] == '\0') break;
  }
  if(_bcondSetName != 0) {
    cinfo->bcondSetName = new char [32];
    for(int i=0; i<32; ++i) {
      if(i==31) cinfo->bcondSetName[i] = '\0';
      else cinfo->bcondSetName[i] = _bcondSetName[i];
       if(_bcondSetName[i] == '\0') break;
    }
  }
}

//---------------------------------------------------------------------

void GeoSource::createSingleCpuToSub(int numSub)
{
  MatrixTimers &mt = domain->getTimers();
  mt.memoryCPUMAP -= memoryUsed();
  int *ptr = new int[2];
  int *trg = new int[numSub];
  ptr[0] = 0;
  ptr[1] = numSub;
  for (int iSub = 0; iSub < numSub; ++iSub)
    trg[iSub] = iSub;
  cpuToSub = new Connectivity(1, ptr, trg);
  mt.memoryCPUMAP += memoryUsed();
}

//---------------------------------------------------------------------

void GeoSource::setControlFile(char *_filename)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->fileName = _filename;
}

//--------------------------------------------------------------------

void GeoSource::setControlRoutine(char *_routinename)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->routineName = _routinename;
}

//--------------------------------------------------------------------

int GeoSource::setSensorLocations(int _numSensor, BCond *_sensor)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numSensor = _numSensor;
  claw->sensor    = _sensor;
  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setActuatorLocations(int _numActuator, BCond *_actuator)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numActuator = _numActuator;
  claw->actuator    = _actuator;
  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setUsddLocation(int _numSensor, BCond *_sensor)
{
  if (claw == 0) claw = new ControlLawInfo;

  claw->numUserDisp = _numSensor;

  // need to copy the pointer data
  claw->userDisp    = new BCond[claw->numUserDisp];
  for (int k = 0; k < claw->numUserDisp; k++)  {
    claw->userDisp[k].nnum = _sensor[k].nnum;
    claw->userDisp[k].dofnum = _sensor[k].dofnum;
    claw->userDisp[k].val = _sensor[k].val;
  }

  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setUsdfLocation(int _numActuator, BCond *_actuator)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numUserForce = _numActuator;
  claw->userForce    = _actuator;
  return 0;
}

//--------------------------------------------------------------------

void GeoSource::outputEnergies(int fileNum, double time, double Wext, double Waero,
			       double Wela, double Wkin, double Wdmp,
                               double error)
{
  fprintf(oinfo[fileNum].filptr," %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
          time, Wext, Waero, Wela, Wkin, Wdmp, error);

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputEnergies(int fileNum, double time, DComplex Wext, DComplex Waero,
                               DComplex Wela, DComplex Wkin, DComplex Wdmp,
                               DComplex error)
{
  // print real part
  fprintf(oinfo[fileNum].filptr," %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
          time, Wext.real(), Waero.real(), Wela.real(), Wkin.real(), Wdmp.real(), error.real());

  // print imaginary part
  fprintf(oinfo[fileNum].filptr," %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
          time, Wext.imag(), Waero.imag(), Wela.imag(), Wkin.imag(), Wdmp.imag(), error.imag());

  fflush(oinfo[fileNum].filptr);
}

//--------------------------------------------------------------------

int GeoSource::getHeaderDescription(char *headDescrip, int fileNumber)
{
  int dataType = 0;  // 1 for nodal, 2 for elemental

  // solver information structure
  SolverInfo sinfo = domain->solInfo();

  char prbType[20];
  if(sinfo.probType == SolverInfo::Static) strcpy(prbType,"Static");
  else if(sinfo.probType == SolverInfo::Dynamic) strcpy(prbType,"Dynam");
  else if(sinfo.probType == SolverInfo::NonLinStatic ||
          sinfo.probType == SolverInfo::MatNonLinStatic) strcpy(prbType,"NLStatic");
  else if(sinfo.probType == SolverInfo::NonLinDynam ||
          sinfo.probType == SolverInfo::MatNonLinDynam) strcpy(prbType,"NLDynamic");
  else if(sinfo.probType == SolverInfo::ArcLength) strcpy(prbType,"Arclength");
  else if(sinfo.probType == SolverInfo::TempDynamic) strcpy(prbType,"Temp");
  else if(sinfo.probType == SolverInfo::AxiHelm) strcpy(prbType,"AxiHelm");
  else if(sinfo.probType == SolverInfo::Modal) strcpy(prbType,"Modal");
  if(isShifted() && sinfo.probType != SolverInfo::Modal) {
    if(sinfo.doFreqSweep) strcpy(prbType,"FrequencySweep");
    else strcpy(prbType,"FrequencyResponse");
  }
  if(sinfo.isAcousticHelm()) {
    if(sinfo.doFreqSweep) strcpy(prbType,"FAcousticSweep");
    else strcpy(prbType,"FAcoustic");
  }

  if (oinfo[fileNumber].type == OutputInfo::YModulus)
    strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::MDensity)
    strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::Thicknes)
    strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::ShapeAtt)
    strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::ShapeStc)
    strcpy(prbType,"Static");

  int type = oinfo[fileNumber].type;

  // No header for CompositeData
  if (type == OutputInfo::Composit)
    return dataType;

  if (type == OutputInfo::InXForce || type == OutputInfo::InYForce ||
      type == OutputInfo::InZForce || type == OutputInfo::AXMoment ||
      type == OutputInfo::AYMoment || type == OutputInfo::AZMoment )  {
    sprintf(headDescrip, header[type], prbType, cinfo->elemSetName, nElem);
    dataType = 2;
    return dataType;
  }

  int avgnum = oinfo[fileNumber].averageFlg;
  if (avgnum == 1 || avgnum == 2)  {
    sprintf(headDescrip, header[type], prbType, cinfo->nodeSetName, numNodes); // PJSA 4-13-05
    dataType = 1;
  }
  else {
    sprintf(headDescrip, ele_header[type], prbType, cinfo->elemSetName);  // PJSA 3-24-05
    dataType = 2;
  }

  return dataType;
}

ControlInterface* GeoSource::getUserSuppliedFunction()
{
  ControlInterface *userSupFunc = 0;

#if !defined(WINDOWS) && !defined(SALINAS)
  if(claw) {
    void *handle;
    dlerror(); // forget about the last error
    handle = dlopen(claw->fileName, RTLD_NOW);
    char *errorMsg;
    if ((errorMsg = dlerror() ) != 0) {
      fprintf(stderr," *** ERROR: in dynamic loading of %s: %s\n",
              claw->fileName,errorMsg);
      exit(-1);
    }

  ControlInterface ** fcp =
      (ControlInterface **) dlsym(handle, claw->routineName);

  if (fcp == 0) {
    fprintf(stderr," *** ERROR: in dynamic loading of %s: "
                     "control function not found\n",
                     claw->routineName);
    exit(-1);
  }

  userSupFunc = *fcp;
  }
#else
	std::cerr << " (W) : Warning : not supported under windows !" << std::endl;
#endif
  return userSupFunc;
}

//--------------------------------------------------------------------

void GeoSource::outputElemStress(int fileNum, double *stressData,
                                 int numOutElems, int *offsets, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time >= 0.0) {
    filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
    if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
  }

  for(int i = 0; i < numOutElems; i++)  {
    int numNodes = offsets[i+1] - offsets[i];
    for(int iNode = 0; iNode < numNodes; iNode++)
      filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode]);
    filePrint(oinfo[fileNum].filptr,"\n");
  }

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputElemStress(int fileNum, DComplex *stressData,
                                 int numOutElems, int *offsets, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  // print real part
  if(time >= 0.0)
    filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, freq());
  for(int i = 0; i < numOutElems; i++)  {
    int numNodes = offsets[i+1] - offsets[i];
    for (int iNode = 0; iNode < numNodes; iNode++)
      filePrint(oinfo[fileNum].filptr," % *.*e", w, p, stressData[offsets[i]+iNode].real());
    filePrint(oinfo[fileNum].filptr,"\n");
  }

  switch(oinfo[fileNum].complexouttype) {
    default:
    case OutputInfo::realimag :
      // print real part
      if(time >= 0.0) {
        filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
      }
      for(int i = 0; i < numOutElems; i++)  {
        int numNodes = offsets[i+1] - offsets[i];
        for(int iNode = 0; iNode < numNodes; iNode++)
          filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode].real());
        filePrint(oinfo[fileNum].filptr,"\n");
      }
      // print imaginary part
      if(time >= 0.0) {
        filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
      }
      for(int i = 0; i < numOutElems; i++)  {
        int numNodes = offsets[i+1] - offsets[i];
        for(int iNode = 0; iNode < numNodes; iNode++)
          filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode].imag());
        filePrint(oinfo[fileNum].filptr,"\n");
      }
      break;
    case OutputInfo::modulusphase :
      // print modulus
      if(time >= 0.0) {
        filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
      }
      for(int i = 0; i < numOutElems; i++)  {
        int numNodes = offsets[i+1] - offsets[i];
        for(int iNode = 0; iNode < numNodes; iNode++)
          filePrint(oinfo[fileNum].filptr," % *.*E", w, p, std::abs(stressData[offsets[i]+iNode]));    //CRW
        filePrint(oinfo[fileNum].filptr,"\n");
      }
      // print phase part
      if(time >= 0.0) {
        filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,time);
        if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
      }
      for(int i = 0; i < numOutElems; i++)  {
        int numNodes = offsets[i+1] - offsets[i];
        for(int iNode = 0; iNode < numNodes; iNode++)
          filePrint(oinfo[fileNum].filptr," % *.*E", w, p, arg(stressData[offsets[i]+iNode]));
        filePrint(oinfo[fileNum].filptr,"\n");
      }
      break;
    case OutputInfo::animate :
      double phi = 0;
      double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
      for(int i=0; i<oinfo[fileNum].ncomplexout; ++i) {
        filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,phi);
        if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
        for(int j = 0; j < numOutElems; j++) {
          int numNodes = offsets[j+1] - offsets[j];
          for(int iNode = 0; iNode < numNodes; iNode++) {
            double proj = std::abs(stressData[offsets[j]+iNode])*cos(arg(stressData[offsets[j]+iNode])-phi);    //CRW
            filePrint(oinfo[fileNum].filptr," % *.*E", w, p, proj);
          }
          filePrint(oinfo[fileNum].filptr,"\n");
        }
        phi += incr;
      }
      break;
  }

  fflush(oinfo[fileNum].filptr);
}

//----------------------------------------------------------------------

int GeoSource::glToPackElem(int i)
{
  if(glToPckElems.find(i) != glToPckElems.end())
    return glToPckElems[i];
  else
    return -1;
}

//-----------------------------------------------------------------------

void
GeoSource::loadMaterial(const char *matName, const char *fileName)
{
#if !defined(WINDOWS) && !defined(SALINAS)
  void *handle;
  handle = dlopen(fileName, RTLD_NOW);
  char *errorMsg;
  if ((errorMsg = dlerror() ) != 0) {
    fprintf(stderr," *** ERROR: in dynamic loading of Material file %s: %s\n",
  	    fileName,errorMsg);
    exit(-1);
  }
  MatLoader fct = (MatLoader) dlsym(handle, "materialLoader");
  if(fct == 0) {
    fprintf(stderr," *** ERROR: in dynamic loading of %s: "
    "materialLoader not found\n", fileName);
    exit(-1);
  }
  fprintf(stderr, "Loaded material is at %p\n", fct);
  char *mm = strdup(matName);
  userDefMat[mm] = fct;

#else
        std::cerr << " (W) : Warning : not supported under windows/salinas !" << std::endl;
#endif
}

void
GeoSource::addMaterial(int i, const char *matName, DoubleList &args)
{
 std::map<const char *, MatLoader, ltstr>::iterator it =
    userDefMat.find(matName);
 if(it == userDefMat.end()) {
    fprintf(stderr, "User defined material %s is not found\n", matName);
    exit(-1);
 }
 MatLoader ml = userDefMat[matName];
 materials[i] = (*ml)(args.nval, args.v);
}

void
GeoSource::simpleDecomposition(int numSubdomains, bool estFlag, bool weightOutFlag)
{
 decJustCalled=true;
 //long m1 = memoryUsed();
 double    t1 = getTime();

 int maxEle = elemSet.last();

 Elemset baseSet(maxEle);
 int nSpring = 0, maxSprNodes = 0, nMass = 0;
 int iEle, iSub;
  for(iEle = 0; iEle < maxEle; ++iEle)
   if(elemSet[iEle]) {
     if(elemSet[iEle]->isSpring() == false && elemSet[iEle]->isMass() == false
        && elemSet[iEle]->isRigidMpcElement() == false)  // PJSA 7-21-05  these are converted to LMPCs for feti solver
       baseSet.elemadd(iEle, elemSet[iEle]);
     else {
       if(elemSet[iEle]->isSpring()) {
	 nSpring++;
	 if(elemSet[iEle]->numNodes() > maxSprNodes)
	   maxSprNodes = elemSet[iEle]->numNodes();
       }
       else nMass++;
     }
   }
 if(nSpring > 0) filePrint(stderr, " ... This mesh has %d springs ...\n", nSpring);

/* already done in GeoSource::addFsiElements
 if ( (domain->solInfo().isCoupled) && !(domain->solInfo().isMatching) &&
      (domain->solInfo().fetiInfo.fsi_corner != 0) ) {
   ResizeArray<LMPCons *> &fsi = domain->getFSI();
// JLchange: see also GeoSource::addFsiElements
//   for(int i = 0; i < domain->getNumFSI(); ++i)
//     baseSet.fsielemadd(maxEle+i, fsi[i]); // PJSA 7-31-06
   int index = maxEle;
   for(int i = 0; i < domain->getNumFSI(); ++i) {
     int fluidNode = fsi[i]->lmpcnum;
     for(int j=0; j< (fsi[i])->nterms; j++) {
       LMPCons *thisFsi = new LMPCons(fluidNode, 0.0);
       LMPCTerm thisLmpcTerm((fsi[i])->terms[j], 1.0);
       thisFsi->addterm(&(thisLmpcTerm));
       baseSet.fsielemadd(index, thisFsi);
       index++;
     }
   }
 }
*/

 MultiFront mf(&baseSet, &nodes, bool(domain->getNumFSI()));
 // filePrint(stderr," ... Constructed Connectivities In %14.5f sec and %14.3f Mb\n", (getTime() - t1)/1000.0,(memoryUsed() - m1)/(1024.0*1024.0));

 // Decompose and optimize the structure into subdomains
 if ( (domain->solInfo().isCoupled) && !(domain->solInfo().isMatching) &&
      (domain->solInfo().fetiInfo.fsi_corner != 0) )
   optDec = mf.decompose(numSubdomains, bool(domain->getNumFSI()));
 else
   optDec = mf.decompose(numSubdomains);

 filePrint(stderr, " ... %d Elements Have Been Arranged in %d Subdomains and %d Springs ...\n",
           optDec->pele[optDec->nsub], optDec->nsub, nSpring);

 nSpring += nMass;
 if(nSpring || nMass) {
   int (*springAssign)[2] = new int[nSpring][2];
   int *subIncr = new int[optDec->nsub];
   for(iSub = 0; iSub < optDec->nsub; ++iSub)
     subIncr[iSub] = 0;
   int *lNd = (int *) dbg_alloca(sizeof(int)*maxSprNodes);
   nSpring = 0;
   for(iEle = 0; iEle < maxEle; ++iEle)
     if(elemSet[iEle] && (elemSet[iEle]->isSpring() || elemSet[iEle]->isMass())) {
       elemSet[iEle]->nodes(lNd);
       springAssign[nSpring][0] = iEle;
       springAssign[nSpring][1] = mf.bestSubFor(elemSet[iEle]->numNodes(), lNd);
       if(springAssign[nSpring][1] < 0)
	 filePrint(stderr, " *** WARNING: Spring element %d is not connected to any non-spring element\n",
                   springAssign[nSpring][0]+1);
       else
	 subIncr[springAssign[nSpring][1]]++;
       nSpring++;
     }
   int *nDecPtr = new int[optDec->nsub+1];
   int *nDecTarget = new int[optDec->pele[optDec->nsub] + nSpring];
   int ptr = 0;
   for(iSub = 0; iSub < optDec->nsub; ++iSub) {
     ptr += subIncr[iSub] + optDec->num(iSub);
     nDecPtr[iSub] = ptr;
   }
   nDecPtr[optDec->nsub] = ptr;
   for(iSub = 0; iSub < optDec->nsub; ++iSub)
     for(iEle = 0; iEle < optDec->num(iSub); ++iEle)
       nDecTarget[--nDecPtr[iSub]] = (*optDec)[iSub][iEle];
   for(iEle = 0; iEle < nSpring; ++iEle)
     nDecTarget[--nDecPtr[springAssign[iEle][1]]] = springAssign[iEle][0];
   delete [] optDec->pele;
   delete [] optDec->eln;
   optDec->pele = nDecPtr;
   optDec->eln  = nDecTarget;

   // fprintf(stderr, "Going to make a check on springs, max = %d\n", maxSprNodes);
   bool *ndIsUsed = new bool[nodes.size()];
   int iNode;
   for(iNode = 0; iNode < nodes.size(); ++iNode)
     ndIsUsed[iNode] = false;
   for(iSub = 0; iSub < optDec->nsub; ++iSub) {
     for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
       int nds[128];
       int enm = (*optDec)[iSub][iEle];
       if(elemSet[enm]->isSpring() == false && elemSet[enm]->isMass() == false) {
	 int nn = elemSet[enm]->numNodes();
	 elemSet[enm]->nodes(nds);
	 for(int i = 0; i < nn; ++i)
	   ndIsUsed[nds[i]] = true;
       }
     }
     for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
       int nds[128];
       int enm = (*optDec)[iSub][iEle];
       if(elemSet[enm]->isSpring()|| elemSet[enm]->isMass()) {
	 elemSet[enm]->nodes(nds);
	 if(ndIsUsed[nds[0]] == false && ndIsUsed[nds[1]] == false)
	   filePrint(stderr, " *** WARNING: Found a badly assigned spring\n");
       }
     }
     for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
       int nds[128];
       int enm = (*optDec)[iSub][iEle];
       if(elemSet[enm]->isSpring() == false && elemSet[enm]->isMass() == false) {
	 int nn = elemSet[enm]->numNodes();
	 elemSet[enm]->nodes(nds);
	 for(int i = 0; i < nn; ++i)
	   ndIsUsed[nds[i]] = false;
       }
     }
   }
   // filePrint(stderr," ... Assigned Springs to Subdomains In %14.5f sec ...\n", (getTime() - t1)/1000.0);
 }

 if(weightOutFlag){
   FILE *weightFile;
   weightFile = domain->openFile(cinfo->checkfile, ".load");
   filePrint(stderr," ... Evaluating Subdomain Weight Distribution ...\n");
   t1 = getTime();
   double *subW = new double[optDec->nsub];
   double minW, maxW = 0;
   double totW = 0;
   int iij;
   for(iij = 0; iij < optDec->nsub; ++iij) {
     double trueW = 0;
     for(iEle = 0; iEle < optDec->num(iij); ++iEle)
       trueW += elemSet[(*optDec)[iij][iEle]]->trueWeight();
     subW[iij] = trueW;
     if(iij == 0 || trueW < minW)
       minW = trueW;
     if(trueW > maxW)
       maxW = trueW;
     totW += trueW;
   }
   fprintf(weightFile, "# %d %e %e %e %f\n", optDec->nsub, minW, totW/optDec->nsub,
   maxW, maxW*optDec->nsub/totW);
   for(iij = 0; iij < optDec->nsub; ++iij)
     fprintf(weightFile, "%e\n", subW[iij]);
   fclose(weightFile);
 }

 if(estFlag) {
   double mem[5];
   // Open and output memory file
   filePrint(stderr," ... Estimating Subdomain Memory Requirements ...\n");
   t1 = getTime();
   FILE *memFilePtr = domain->openFile(cinfo->checkfile, ".mem");
   mf.memEstimate(optDec, 3, mem, memFilePtr);
   fclose(memFilePtr);
 }

 // Open optimized decomposition file
 filePrint(stderr," ... Saving Decomposition File      ...\n");
 t1 = getTime();
 FILE *optFilePtr = domain->openFile(cinfo->checkfile, ".optDec");

 // Output optimized decomposition to file
 optDec->setName(cinfo->elemSetName);
#ifndef SALINAS
 optDec->outputDump(optFilePtr, 0);
#endif
 fclose(optFilePtr);
}

void
GeoSource::setExitAfterDec(bool exit)
{
  exitAfterDec=exit;
  if(exit) domain->solInfo().setProbType(SolverInfo::Decomp);
}

#ifndef SALINAS
#include <Driver.d/Sower.h>
#include <vector>
/*
 * writeDistributedInputFiles will do the first step of the distributed input for fem :
 * it will assign subdomains to nCluster clusters and create a data file for each cluster
 * containing only the data pertinent to that cluster. (all of this is done by creating and initializing a Sower object properly )
 * @param nCluster number of cluster fem will run on
 * @see class Sower Driver.d/Sower.h
*/
//#define SOWER_DEBUG
void GeoSource::writeDistributedInputFiles(int nCluster, Domain *domain)
{
  if(subToElem == 0)
    getTextDecomp(true); //the option true or false is for sowering or not
                         // if sowering there is no reordering

  // TOPOLOGY (also element data PRESSURE and PRELOAD)
  Sower sower(subToElem, domain->getElementSet(), nCluster, domain->viewSurfEntities()); //HB
  /* the sower object creates the cluster to element connectivity and from it will
     create cluster to data connectivities for each data added subsequently.
     Thus It will finally know which cluster needs what */

  // begin adding data to sower object

  // NODES
  Connectivity *eToN = new Connectivity(&(domain->getElementSet()));
  sower.addParentToChildData<NodesIO, CoordSet*, Connectivity*>(NODES_TYPE, ELEMENTS_TYPE, 0, &nodes, eToN);
  delete eToN;

  // ATTRIBUTES
  std::pair<int, map<int, Attrib>* > attrPair = std::make_pair(na, &attrib);
  ImplicitConnectivity<std::pair<int,map<int, Attrib>* >*, ElemAttrAccessor>
    *EtoAtt = new ImplicitConnectivity<std::pair<int, map<int, Attrib>* >*, ElemAttrAccessor>(&attrPair);
  sower.addParentToChildData<AttribIO>(ATTRIBUTES_TYPE, ELEMENTS_TYPE, 0, &attrPair, EtoAtt);

  // NEUMAN BOUNDARY CONDITIONS
  typedef std::pair<int, BCond*> BCPair;
  typedef ImplicitConnectivity<BCPair*, BCDataAccessor> implicitBC;
  BCPair forcesPair = std::make_pair(numNeuman, nbc);
  implicitBC *forcesToNodes = new implicitBC(&forcesPair);
  sower.addChildToParentData<BCDataIO<FORCES_TYPE> >(FORCES_TYPE, NODES_TYPE, 0, &forcesPair, forcesToNodes);
  delete forcesToNodes;

  // DIRICHLET BOUNDARY CONDITIONS
  BCPair dispPair = std::make_pair(numDirichlet, dbc);
  implicitBC *dispToNodes = new implicitBC(&dispPair);
  sower.addChildToParentData<BCDataIO<DISPLACEMENTS_TYPE> >(DISPLACEMENTS_TYPE, NODES_TYPE, 0, &dispPair, dispToNodes);
  delete dispToNodes;

  // MATERIALS
  typedef std::pair<int, SPropContainer* > MATPair;
  typedef ImplicitConnectivity<std::pair<int,map<int,Attrib>* >*, MatAttrAccessor> implicitMat;
  MATPair matPair = std::make_pair(numProps, &sProps);
  implicitMat* EleToMat = new implicitMat(&attrPair);
  sower.addParentToChildData<MatIO>(MATERIALS_TYPE, ELEMENTS_TYPE, 0, &matPair, EleToMat);
  delete EleToMat;

  delete EtoAtt; // not used any further

  // TETT
  typedef ImplicitConnectivity<std::pair<int,SPropContainer* >*, CurveMatAccessor> implicitCurve;
  implicitCurve* MatToCurve = new implicitCurve(&matPair);
  sower.addParentToChildData<MFTTDataIO<TETT_TYPE> >(TETT_TYPE, MATERIALS_TYPE, 0, domain->getCTETT(), MatToCurve);
  delete MatToCurve;

  // YMTT
  typedef ImplicitConnectivity<std::pair<int,SPropContainer* >*, CurveYoungMatAccessor> implicitYoungCurve;
  implicitYoungCurve* MatYoungToCurve = new implicitYoungCurve(&matPair);
  sower.addParentToChildData<MFTTDataIO<YMTT_TYPE> >(YMTT_TYPE, MATERIALS_TYPE, 0, domain->getYMTT(), MatYoungToCurve);
  delete MatYoungToCurve;

  // LMPC
  typedef ImplicitConnectivity<std::pair<int,ResizeArray<LMPCons *>* >*, LMPCAccessor> implicitLMPC;
  std::pair<int,ResizeArray<LMPCons *>* > LMPCPair = std::make_pair(domain->getNumLMPC(),domain->getLMPC());
  implicitLMPC* ImplicitLMPC = new implicitLMPC(&LMPCPair);
  sower.addChildToParentData<LMPCIO>(LMPC_TYPE, NODES_TYPE, 0, &LMPCPair, ImplicitLMPC);
  delete ImplicitLMPC;

  // COMPOSITE
  typedef std::pair<int,map<int,Attrib>* > AttPair;
  typedef ImplicitConnectivity<AttPair*, CmpAttrAccessor> implicitCMP;
  implicitCMP* EleToCmp =  new implicitCMP(&attrPair);
  typedef std::pair<int, ResizeArray<LayInfo *>* > layInfoPair;
  layInfoPair  lip = std::make_pair(numLayInfo, &layInfo);
  sower.addParentToChildData<CompositeLIO>(COMPOSITEL_TYPE, ELEMENTS_TYPE, 0, &lip, EleToCmp);
  typedef std::pair<int, ResizeArray<CoefData *>* > coefDataPair;
  coefDataPair cdp = std::make_pair(numCoefData, &coefData);
  sower.addParentToChildData<CompositeCIO>(COMPOSITEC_TYPE, ELEMENTS_TYPE, 0, &cdp, EleToCmp);

  // CFRAMES
  typedef ImplicitConnectivity<AttPair*, CmpFrAttrAccessor> implicitCFM;
  implicitCFM* EleToCFM =  new implicitCFM(&attrPair);
  typedef std::pair<int, ResizeArray<double *>* > cFramesPair;
  cFramesPair cfp = std::make_pair(numCframes, &cframes);
  sower.addParentToChildData<CFramesIO>(CFRAMES_TYPE, ELEMENTS_TYPE, 0, &cfp, EleToCFM);

  // BOFFSET
  typedef ImplicitConnectivity<vector<OffsetData>*, BoffsetAccessor> implicitBoffset;
  implicitBoffset* offToElem = new implicitBoffset(&offsets);
  sower.addChildToParentData<BoffsetIO>(BOFFSET_TYPE, ELEMENTS_TYPE, 0, &offsets, offToElem);
  delete offToElem;

  // EFRAMES
  typedef std::pair<int,ResizeArray<EFrameData>* > EfPair;
  typedef ImplicitConnectivity<EfPair*, EFrameDataAccessor > implicitEFrame;
  EfPair efPair = std::make_pair(numEframes,&efd);
  implicitEFrame* efToElem = new implicitEFrame(&efPair);
  sower.addChildToParentData<EFrameIO>(EFRAME_TYPE, ELEMENTS_TYPE, 0, &efPair, efToElem);

  // DIMASS
  vector<DMassData*> dmv; // what can I do with a chained list !
  DMassData* fdm = domain->getFirstDMassData();
  if(fdm != 0) // some DMASS
    {
      do
	{
	  dmv.push_back(fdm);
	  fdm = fdm->next;
	}
      while(fdm->next!= 0);
      typedef ImplicitConnectivity<vector<DMassData* >*, DimassAccessor > implicitDMass;
      implicitDMass * massToNode = new implicitDMass(&dmv);
      sower.addChildToParentData<DMassIO>(DIMASS_TYPE, NODES_TYPE, 0, &dmv, massToNode);
    }

  // CONVECTION
  BCPair convPair = std::make_pair(numConvBC, cvbc);
  implicitBC *convToNodes = new implicitBC(&convPair);
  sower.addChildToParentData<BCDataIO<CONV_TYPE> >(CONV_TYPE, NODES_TYPE, 0, &convPair, convToNodes);
  delete convToNodes;

  // RADIATION
  BCPair radPair = std::make_pair(numRadBC, rdbc);
  implicitBC *radToNodes = new implicitBC(&radPair);
  sower.addChildToParentData<BCDataIO<RAD_TYPE> >(RAD_TYPE, NODES_TYPE, 0, &radPair, radToNodes);
  delete radToNodes;

  // INITIAL DISPLACEMENTS
  BCPair idispPair = std::make_pair(numIDis, iDis);
  implicitBC *idispToNodes = new implicitBC(&idispPair);
  sower.addChildToParentData<BCDataIO<IDISP_TYPE> >(IDISP_TYPE, NODES_TYPE, 0, &idispPair, idispToNodes);
  delete idispToNodes;

  // INITIAL DISPLACEMENTS (6 column)
  BCPair idisp6Pair = std::make_pair(numIDis6, iDis6);
  implicitBC *idisp6ToNodes = new implicitBC(&idisp6Pair);
  sower.addChildToParentData<BCDataIO<IDISP6_TYPE> >(IDISP6_TYPE, NODES_TYPE, 0, &idisp6Pair, idisp6ToNodes);
  delete idisp6ToNodes;

  // INITIAL VELOCITIES
  BCPair ivelPair = std::make_pair(numIVel, iVel);
  implicitBC *ivelToNodes = new implicitBC(&ivelPair);
  sower.addChildToParentData<BCDataIO<IVEL_TYPE> >(IVEL_TYPE, NODES_TYPE, 0, &ivelPair, ivelToNodes);
  delete ivelToNodes;

  // INITIAL TEMPERATURES
  BCPair itempPair = std::make_pair(numITemp, iTemp);
  implicitBC *itempToNodes = new implicitBC(&itempPair);
  sower.addChildToParentData<BCDataIO<ITEMP_TYPE> >(ITEMP_TYPE, NODES_TYPE, 0, &itempPair, itempToNodes);
  delete itempToNodes;

  // COMPLEX DIRICHLET
  typedef std::pair<int, ComplexBCond*> ComplexBCPair;
  typedef ImplicitConnectivity<ComplexBCPair*, ComplexBCDataAccessor> implicitComplexBC;
  ComplexBCPair cdirPair = std::make_pair(numComplexDirichlet, cdbc);
  implicitComplexBC *cdirToNodes = new implicitComplexBC(&cdirPair);
  sower.addChildToParentData<ComplexBCDataIO<HDIR_TYPE> >(HDIR_TYPE, NODES_TYPE, 0, &cdirPair, cdirToNodes);
  delete cdirToNodes;

  // COMPLEX NEUMANN
  ComplexBCPair cneuPair = std::make_pair(numComplexNeuman, cnbc);
  implicitComplexBC *cneuToNodes = new implicitComplexBC(&cneuPair);
  sower.addChildToParentData<ComplexBCDataIO<HNEU_TYPE> >(HNEU_TYPE, NODES_TYPE, 0, &cneuPair, cneuToNodes);
  delete cneuToNodes;

  // ATDDNB & HDNB
  typedef std::pair<int, SommerElement **> SommerPair;
  SommerPair dnbPair = std::make_pair(domain->numNeum,(domain->neum).yield());
  Connectivity *dnbToE = new Connectivity( &domain->getElementSet(), dnbPair.first, dnbPair.second);
  sower.addChildToParentData<SommerDataIO<DNB_TYPE> >(DNB_TYPE, ELEMENTS_TYPE, 0, &dnbPair, dnbToE);
  delete dnbToE;

  // ATDROB & HSCB
  SommerPair scatPair = std::make_pair(domain->numScatter,(domain->scatter).yield());
  Connectivity *scatToE = new Connectivity(&domain->getElementSet(), scatPair.first, scatPair.second);
  sower.addChildToParentData<SommerDataIO<SCAT_TYPE> >(SCAT_TYPE, ELEMENTS_TYPE, 0, &scatPair, scatToE);
  delete scatToE;

  // ATDARB & HARB
  SommerPair arbPair = std::make_pair(domain->numSommer,(domain->sommer).yield());
  Connectivity *arbToE = new Connectivity(&domain->getElementSet(), arbPair.first, arbPair.second);
  sower.addChildToParentData<SommerDataIO<ARB_TYPE> >(ARB_TYPE, ELEMENTS_TYPE, 0, &arbPair, arbToE);
  delete arbToE;

  if (claw) {
    // distribute control law data -> have to do that for all 4 categories...
    for (int i = 0 ; i < claw->numSensor ; i++) {
      claw->sensor[i].val=i;
    }
    BCPair sensorPair = std::make_pair(claw->numSensor, claw->sensor);
    implicitBC *sensorToNodes = new implicitBC(&sensorPair);
    sower.addChildToParentData<BCDataIO<SENSOR_TYPE> >(SENSOR_TYPE, NODES_TYPE, 0, &sensorPair, sensorToNodes);
    delete sensorToNodes;

    for (int i = 0 ; i < claw->numActuator ; i++) {
      claw->actuator[i].val=i;
    }
    BCPair actuatorPair = std::make_pair(claw->numActuator, claw->actuator);
    implicitBC *actuatorToNodes = new implicitBC(&actuatorPair);
    sower.addChildToParentData<BCDataIO<ACTUATOR_TYPE> >(ACTUATOR_TYPE, NODES_TYPE, 0, &actuatorPair, actuatorToNodes);
    delete actuatorToNodes;

    for (int i = 0 ; i < claw->numUserDisp ; i++) {
      claw->userDisp[i].val=i;
    }
    BCPair usddPair = std::make_pair(claw->numUserDisp, claw->userDisp);
    implicitBC *usddToNodes = new implicitBC(&usddPair);
    sower.addChildToParentData<BCDataIO<USDD_TYPE> >(USDD_TYPE, NODES_TYPE, 0, &usddPair, usddToNodes);
    delete usddToNodes;

    for (int i = 0 ; i < claw->numUserForce ; i++) {
      claw->userForce[i].val=i;
    }
    BCPair usdfPair = std::make_pair(claw->numUserForce, claw->userForce);
    implicitBC *usdfToNodes = new implicitBC(&usdfPair);
    sower.addChildToParentData<BCDataIO<USDF_TYPE> >(USDF_TYPE, NODES_TYPE, 0, &usdfPair, usdfToNodes);
    delete usdfToNodes;
  }

  // done adding data to sower object
#ifdef SOWER_DEBUG
  cerr << "  Debug requested for sower\n" << endl;
  sower.printDebug();
#endif
  sower.write();

}
#endif

Connectivity *
GeoSource::getDecomposition()
{
  if(binaryInput) getBinaryDecomp();
  else if(!optDec) getTextDecomp();
  else {
    int i;
    int numSub = optDec->nsub;
    // Allocate memory for offsets
    int *cx = new int[numSub+1];
    // Allocate memory for subdomain to element connectivity
    nElem=0;
    for(i=0; i<numSub; i++) nElem += optDec->num(i);
    int *connect = new int[nElem];
    for(i=0; i<=numSub; i++) cx[i] = optDec->pele[i];
    for(i=0; i<nElem; i++) connect[i] = optDec->eln[i];
    subToElem = new Connectivity(numSub,cx,connect);
    subToElem->renumberTargets(glToPckElems);  // PJSA: required if gaps in element numbering

#ifdef DISTRIBUTED // PJSA 1-22-07
    int* ptr = new int[numSub+1];
    int* target = new int[numSub];
    ptr[0] = 0;
    for(int i=0;i<numSub; i++){
       ptr[i+1] = ptr[i]+1;
       target[i] = 0;
    }
    subToClus = new Connectivity(numSub,ptr,target);
    numClusters = 1;
    clusToSub = subToClus->reverse();
    numClusNodes = nGlobNodes;
    numClusElems = nElem; //HB: not sure this is be always correct (i.e. phantoms els ...)
#endif


  }
  return subToElem;
}

void GeoSource::getBinaryDecomp()
{
  if(!subToElem) {
    BinFileHandler fp("decomposition", "rb");
    subToElem = new Connectivity(fp);
    subToNode = new Connectivity(fp);
#ifdef SOWER_DEBUG
    cerr << "*** subToElem, from decomposition binary file: \n"; subToElem->print();
    cerr << "*** subToNode, from decomposition binary file: \n"; subToNode->print();
#endif
  }
}

#ifndef SALINAS
void GeoSource::readGlobalBinaryData()
{
  if(!subToSub || !subToClus) {
    BinFileHandler fp2("connectivity", "rb");
    clusToSub = new Connectivity(fp2);
    numClusters = clusToSub->csize();
    subToClus = clusToSub->reverse();
    subToSub = new Connectivity(fp2);
#ifdef SOWER_DEBUG
    cerr << "*** subToClus, from connectivity binary file: \n"; subToClus->print();
    cerr << "*** subToSub, from connectivity binary file: \n"; subToSub->print();
#endif

    // build global to cluster subdomain map
    int numSub = subToSub->csize();
    gl2ClSubMap = new int[numSub];
    for(int iClus = 0; iClus < numClusters; iClus++)  {
      int clusNum = 0;
      for(int iSub = 0; iSub < clusToSub->num(iClus); iSub++)
        gl2ClSubMap[ (*clusToSub)[iClus][iSub] ] = clusNum++;
    }

    fp2.read(&nGlobNodes, 1);
    domain->setNumNodes(nGlobNodes);
    int nGlobElems;
    fp2.read(&nGlobElems, 1);
    domain->setNumElements(nGlobElems);
#ifdef SOWER_DEBUG
    cerr << "*** global number of nodes = " << nGlobNodes << endl;
    cerr << "*** global number of elements = " << nGlobElems << endl;
#endif
    if(numClusters == 1) {
      numClusNodes = nGlobNodes;
      numClusElems = nGlobElems;
    }
    else {
      filePrint(stdout," *** ERROR: only one cluster is supported \n");
      exit(-1);
    }

    int nGlobSurfs;
    fp2.read(&nGlobSurfs,1);
    domain->setNumSurfs(nGlobSurfs);
#ifdef SOWER_DEBUG
    cerr << " *** global number of surfaces = " << nGlobSurfs << endl;
#endif
    if(nGlobSurfs > 0) {
      Sower sower;
      BinFileHandler *f = sower.openBinaryFile(0);
      for(int isurf=0; isurf<domain->getNumSurfs(); isurf++) {
        int* dummy= 0;
        SurfaceEntity* surf = sower.read<SurfaceIO>(*f, isurf, dummy, true);
        if(dummy) { delete [] dummy; dummy= 0; } // not used
        domain->AddSurfaceEntity(surf,isurf);
      }
    }
  }
}
#endif

void GeoSource::computeClusterInfo(int glSub)
{
  int clusNum = (*subToClus)[glSub][0];
  Connectivity *clusToNode = clusToSub->transcon(subToNode);
  numClusNodes = clusToNode->num(clusNum);
  delete clusToNode;
  Connectivity *clusToElem = clusToSub->transcon(subToElem);
  numClusElems = clusToElem->num(clusNum);
  delete clusToElem;
}

void GeoSource::makeEframe(int ele, int refnode, double *d)
{
  d[0] = d[1] = d[2] = 0.0;
  int beam_nodes[2];
  elemSet[ele]->nodes(beam_nodes);
  Node *node1 = nodes[beam_nodes[0]];
  Node *node2 = nodes[refnode-1];
  d[3] = node2->x - node1->x;
  d[4] = node2->y - node1->y;
  d[5] = node2->z - node1->z;
  d[6] = d[7] = d[8] = 0.0;
}


void GeoSource::setGroupAttribute(int a, int g)
{
 group[g].attributes.push_back(a);
}


void GeoSource::printGroups()
{
  for(map<int, Group >::iterator it = group.begin(); it != group.end(); ++it) {
    cerr << "group_id = " << it->first+1 << endl;
    cerr << "  number of attributes in this group =  " << it->second.attributes.size() << endl;
    cerr << "  number of random properties in this group =  " << it->second.randomProperties.size() << endl;
    cerr << "  attributes: ";
    for(int i = 0; i < int(it->second.attributes.size()); ++i) cerr << it->second.attributes[i]+1 << " ";
    cerr << "\n  random properties: ";
    for(int i = 0; i < int(it->second.randomProperties.size()); ++i) {
      cerr << "rprop = " << it->second.randomProperties[i].rprop
           << ", mean = " << it->second.randomProperties[i].mean
           << ", std_dev = " << it->second.randomProperties[i].std_dev << endl;
     }
    cerr << endl;
  }
}

void GeoSource::setGroupRandomProperty(int g, Rprop prop_type, double mean, double std_dev)
{
#ifndef SALINAS
  RandomProperty rp(prop_type, mean, std_dev);
  group[g].randomProperties.push_back(rp);
  sfem->updatendim();
#endif
}

bool GeoSource::elemOutput()
{
  for(int iInfo = 0; iInfo < geoSource->getNumOutInfo(); iInfo++) {
    if(oinfo[iInfo].averageFlg == 0) return true;
  }
  return false;
}

ControlLawInfo::ControlLawInfo()
{
  fileName = routineName = 0;
  numSensor = numActuator = numUserDisp = numUserForce = 0;
  sensor = 0; actuator = 0; userDisp = 0; userForce = 0;
}

ControlLawInfo::~ControlLawInfo()
{
  if(sensor) delete [] sensor;
  if(actuator) delete [] actuator;
  if(userDisp) delete [] userDisp;
  if(userForce) delete [] userForce;
}

void ControlLawInfo::print()
{
  fprintf(stderr, " Number of Sensors: %d\n", numSensor);
  fprintf(stderr, " Number of Actuators: %d\n", numActuator);
  fprintf(stderr, " Number of UserForces: %d\n", numUserForce);
  fprintf(stderr, " Number of UserDisps: %d\n", numUserDisp);
  fprintf(stderr, " filename: %s\n", fileName);
  fprintf(stderr, " routine: %s\n", routineName);

  for (int j = 0; j < numUserDisp; j++)
    fprintf(stderr, " usdd[%d][%d] = %f\n", userDisp[j].nnum, userDisp[j].dofnum, userDisp[j].val);
}

void ControlLawInfo::makeGlobalClaw(ControlLawInfo *subClaw)
{
  //only the numbers are augmented for sower input
  numSensor += subClaw->numSensor;
  numActuator += subClaw->numActuator;
  numUserForce += subClaw->numUserForce;
  numUserDisp += subClaw->numUserDisp;
}


