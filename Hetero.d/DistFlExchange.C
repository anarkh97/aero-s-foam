#include <Utils.d/dbg_alloca.h>
#include <Utils.d/Connectivity.h>
#include <Math.d/SparseMatrix.h>
#include <Hetero.d/FlExchange.h>
#include <Hetero.d/DistFlExchange.h>
#include <Hetero.d/FilteredFile.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/SubDomain.h>

#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Feti.d/DistrVector.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/SysState.h>

#include <Corotational.d/DistrGeomState.h>

#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

//--------------------------------------------------------------------

// global variables
extern Communicator *structCom, *fluidCom, *heatStructCom;

int toFluid = 1;
int fromFluid = 1;
int toStruc = 2;
int fromTemp = 3;

//---------------------------------------------------------------------

DistFlExchanger::DistFlExchanger(CoordSet **_cs, Elemset **_eset, 
	                         DofSetArray **_cdsa, DofSetArray **_dsa,
                                 OutputInfo *_oinfo)
 : cs(_cs), eset(_eset), cdsa(_cdsa), dsa(_dsa), oinfo(_oinfo),
   tmpDisp(0), dsp(0), vel(0), acc(0), pVel(0)
{
  useFaceElem = false;
}

DistFlExchanger::DistFlExchanger(CoordSet **_cs, Elemset **_eset, SurfaceEntity *_surface, CoordSet *_globalCoords,
                                 Connectivity *_nodeToElem, Connectivity *_elemToSub, SubDomain **_sd,
                                 DofSetArray **_cdsa, DofSetArray **_dsa, OutputInfo *_oinfo)
 : cs(_cs), eset(_eset), surface(_surface), globalCoords(_globalCoords), nodeToElem(_nodeToElem),
   elemToSub(_elemToSub), sd(_sd), cdsa(_cdsa), dsa(_dsa), oinfo(_oinfo),
   tmpDisp(0), dsp(0), vel(0), acc(0), pVel(0)
{
  useFaceElem = true;
}

DistFlExchanger::~DistFlExchanger()
{
  if(useFaceElem) {
    for(int i = 0; i < geoSource->getCpuToSub()->numConnect(); ++i)
      delete fnId2[i];
    delete [] fnId2;
  }
}

//---------------------------------------------------------------------

MatchMap*
DistFlExchanger::getMatchData()
{
  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  MatchMap* globMatches = new MatchMap();

  // create global sub to local sub map
  int totSub = cpuToSub->numConnect();
  int *glToLocSubMap = new int[totSub];
 
  for(int iSub = 0; iSub < totSub; iSub++)
    glToLocSubMap[iSub] = -1;

  for(int iSub = 0; iSub < numSub; iSub++)
    glToLocSubMap[ (*cpuToSub)[myCPU][iSub] ] = iSub;


  if(!useFaceElem) {
    for(int iSub = 0; iSub < numSub; iSub++) {

      int glSub = (*cpuToSub)[myCPU][iSub];
      int locSub = glToLocSubMap[glSub];
      if(locSub < 0)
        fprintf(stderr,"*** ERROR: Bad Global to Local Sub Map\n");

      MatchData *matchData = geoSource->getMatchData(locSub);
      int *numMatch = geoSource->getNumMatchData();

      for(int iMatch = 0; iMatch < numMatch[locSub]; iMatch++) {
        int gPos = matchData[iMatch].glPos;
        InterpPoint iPoint;
        iPoint.subNumber = iSub;
        iPoint.elemNum = matchData[iMatch].elemNum;
        iPoint.xy[0] = matchData[iMatch].xi;
        iPoint.xy[1] = matchData[iMatch].eta;
        iPoint.gap[0] = 0.0;
        iPoint.gap[1] = 0.0;
        iPoint.gap[2] = 0.0;
        (*globMatches)[gPos] = iPoint;
      }
    }
  }
  else {
    // info about surface
    FaceElemSet& feset = surface->GetFaceElemSet();
    CoordSet&   fnodes = surface->GetNodeSet();
    int*          fnId = surface->GetPtrGlNodeIds();

    // assign face elements to a subdomain 
    int totele = elemToSub->csize();
    int *eleTouch = new int[totele];
    int *eleCount = new int[totele];
    for(int j = 0; j < totele; j++) eleTouch[j] = -1;
    int numFaceEle = feset.last();
    int *faceElemToSub = new int[numFaceEle];
    for(int j = 0; j < numFaceEle; ++j) {
      int glEle = feset[j]->findEle(nodeToElem, eleTouch, eleCount, j, fnId);
      int glSub = (*elemToSub)[glEle][0];
      faceElemToSub[j] = glToLocSubMap[glSub];
      
    }
    delete [] eleTouch;
    delete [] eleCount;

    // find node to elem connectivity and local coords. 
    // local coords convention given in "Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.C"
    map<int,locoord> exy = feset.computeNodeLocalCoords(fnId, fnodes.size());
    map<int,locoord>::iterator it;
    for(int j = 0; j < fnodes.size(); ++j) {
      it = exy.find(j);
      if(it==exy.end()) { fprintf(stderr,"Oh no!\n"); exit(-1); }

      int faceElemNum = (it->second).first;
      int subNumber = faceElemToSub[faceElemNum];
      if(subNumber > -1) {
        InterpPoint iPoint;
        iPoint.subNumber = subNumber;
        iPoint.elemNum = faceElemNum;
        iPoint.xy[0] = (it->second).second.first;
        iPoint.xy[1] = (it->second).second.second;
        iPoint.gap[0] = 0.0;
        iPoint.gap[1] = 0.0;
        iPoint.gap[2] = 0.0; /*KW:no gap*/
        (*globMatches)[j] = iPoint;
      }
    }
    delete [] faceElemToSub;

    fnId2 = new int * [numSub];
    for(int i = 0; i < numSub; ++i) {
      fnId2[i] = new int[fnodes.size()];
      for(int j = 0; j < fnodes.size(); ++j) fnId2[i][j] = sd[i]->globalToLocal(fnId[j]);
    }
  }

  delete [] glToLocSubMap;

  return globMatches;
}

//----------------------------------------------------------------------

// This routine negotiate with the fluid codes which match points go where
void DistFlExchanger::negotiate()
{
  FaceElemSet *feset;
  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
    sendEmbeddedWetSurface();
  }
  
  int numFl = fluidCom->remoteSize();  // number of fluid mpi processes
  int *flSize = new int[numFl];  // number of matches per fluid mpi 

  int iFluid;
  for (iFluid = 0; iFluid < numFl; ++iFluid) {
    int tag = FL_NEGOT;
    int nFlMatched;
    int buffLen = 1;
    //int rsize; 
    //int zero = 0;
    RecInfo rInfo = fluidCom->recFrom(tag, &nFlMatched, buffLen);
    int fromNd = rInfo.cpu;
    flSize[fromNd] = nFlMatched;
  }

  // allocate with number of fluid cpu's doing a negotiate
  int *consOrigin = new int[numFl];
  nbSendTo = new int[numFl];

  int nSender = 0;  // number of fluid mpi's sending data

  int maxSize = 0;
  for (iFluid = 0; iFluid < numFl; ++iFluid) {

    if (flSize[iFluid] > 0) {

      nbSendTo[iFluid] = flSize[iFluid];
      consOrigin[iFluid] = nSender;
      if (flSize[iFluid] > maxSize) 
        maxSize = flSize[iFluid];
      nSender++;
      
    }
  }

  MatchMap* globMatches = getMatchData();

  // # of fluid mpi's receiving data from this mpi
  nbrReceivingFromMe = 0;

  int *index = new int[maxSize];
  //fprintf(stderr,"struct %d allocating sndTable of size %d\n", structCom->myID(), nSender);
  sndTable = new InterpPoint *[nSender];

  // receive matches from each sender
  int bufferLen = maxSize;

  int totSize = 0;  // total number of matches received by structure
  int maxNDof = 0;  // to set size of local force vector

  for (iFluid = 0; iFluid < nSender; ++iFluid) {

    int tag = FL_NEGOT+1;
    RecInfo rInfo = fluidCom->recFrom(tag, index, bufferLen);
    int fromNd = rInfo.cpu;
    //int rsize = rInfo.len;
    //int sender = consOrigin[fromNd];  // gives fluid mpi # TODO check if this is correct see FlExchange.C
    int sender = rInfo.cpu;
    
    //fprintf(stderr,"Struct cpu %d receives %d matches from fluid %d\n", structCom->myID(), nbSendTo[sender], sender); 

    // determine how many matches are in this mpi process
    int numHits = 0;
    int iMatch;
    for (iMatch = 0; iMatch < nbSendTo[sender]; iMatch++)
      if (globMatches->find(index[iMatch]) != globMatches->end()) 
        numHits++;

    totSize += numHits;
    int *matchList;  // list of matches

    int totalNDof = 0;

    //fprintf(stderr,"Struct CPU %d has %d matches w/fluid %d\n", structCom->myID(), numHits, sender);
    matchList = new int[numHits];
    sndTable[sender] = new InterpPoint[numHits];

    if (numHits > 0)  {
      // create list of matches in this mpi
      int iHit = 0;
      for (iMatch = 0; iMatch < nbSendTo[sender]; iMatch++)
        if (globMatches->find(index[iMatch]) != globMatches->end()) 
          matchList[iHit++] = iMatch;

      // reset nbSendTo to actual number of matches in this mpi 
      nbSendTo[sender] = numHits;

      // assign match data to the sndTable
      //for (int ipt = 0; ipt < numHits; ++ipt) {
      for (int ipt = 0; ipt < nbSendTo[sender]; ++ipt) {
       
        int arrayPos = index[matchList[ipt]];
        sndTable[sender][ipt].subNumber = (*globMatches)[arrayPos].subNumber; 
        sndTable[sender][ipt].elemNum = (*globMatches)[arrayPos].elemNum; 
        sndTable[sender][ipt].xy[0] = (*globMatches)[arrayPos].xy[0]; 
        sndTable[sender][ipt].xy[1] = (*globMatches)[arrayPos].xy[1]; 
        sndTable[sender][ipt].gap[0] = 0;
        sndTable[sender][ipt].gap[1] = 0;
        sndTable[sender][ipt].gap[2] = 0;

        // assign dofs to element 
        int locSub = sndTable[sender][ipt].subNumber;
        int locElem = sndTable[sender][ipt].elemNum;
        int nDof;

        if(!useFaceElem) {
          Element *thisElement = (*eset[locSub])[locElem];
          if(thisElement == NULL) {
            cerr << " ERROR: DistFlExchanger::negotiate() cpu " << structCom->myID() << ", locSub = " 
                 << locSub << ", locElem = " << locElem << endl;
          }
          nDof = thisElement->numDofs();
        }
        else {
          FaceElement *thisFaceElement = (*feset)[locElem];
          nDof = thisFaceElement->numDofs();
        }
        totalNDof += nDof;
        if (nDof > maxNDof) 
          maxNDof = nDof;
      }
      nbrReceivingFromMe++;
    }
/*
    else  {
      sndTable[sender] = 0;
      matchList = 0;
    }
*/

    // allocate the local dofset array
    int *array = new int[totalNDof];

    // populate the dofs in the send table
    for (int iData = 0; iData < numHits; iData++)  {

        int locSub = sndTable[sender][iData].subNumber;
        int locElem = sndTable[sender][iData].elemNum;

        sndTable[sender][iData].dofs = array;

        if(!useFaceElem) {
          Element *thisElement = (*eset[locSub])[locElem];
          thisElement->dofs(*(cdsa[locSub]), array);
          array += thisElement->numDofs();
        }
        else {
          FaceElement *thisFaceElement = (*feset)[locElem];
          thisFaceElement->dofs(*(cdsa[locSub]), array, fnId2[locSub]); // TODO fnId2 must map from embedded surface to locSub node numbers
          array += thisFaceElement->numDofs();
        }

    }

    tag = FL_NEGOT;
    int one = 1;
    fluidCom->sendTo(fromNd, tag, &numHits, one);

    if (numHits > 0)  {
      tag = FL_NEGOT+1;
      fluidCom->sendTo(fromNd, tag, matchList, numHits);
    }

    // To make sure we can reuse the buffers
    fluidCom->waitForAllReq();

  }

  // allocate the size of local force vector
  localF = new double[maxNDof];

  // create list of fluid mpi's to send to
  idSendTo = new int[nbrReceivingFromMe]; 

  int count = 0;
  for (iFluid = 0; iFluid < nSender; iFluid++)
    if (sndTable[iFluid] != 0)
      idSendTo[count++] = iFluid;

  delete [] index;
  buffer = new double[6*totSize];
  setBufferLength(6*totSize);
}

//---------------------------------------------------------------------

void
DistFlExchanger::sendDisplacements(SysState<DistrVector>& dState, 
                                   double** usrDefDisps, double** usrDefVels, int tag,
                                   DistrGeomState *distrGeomState)
{
  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  if (tmpDisp == 0)
    tmpDisp = new DistrVector(dState.getDisp());

  *tmpDisp = dState.getDisp();
  tmpDisp->linAdd(dt*alpha[0], dState.getVeloc(), dt*alpha[1], dState.getPrevVeloc());
  //if(verboseFlag)
  //  filePrint(stderr, "Disp Norm %e Veloc Norm %e\n", (*tmpDisp)*(*tmpDisp), dState.getVeloc()*dState.getVeloc());

  FaceElemSet *feset;
  int         *fnId;
  if(useFaceElem){
    feset = &(surface->GetFaceElemSet());
    fnId  = surface->GetPtrGlNodeIds();
  }

  Element     *thisElement;
  FaceElement *thisFaceElem;

  double xxx = 0, yyy = 0;
  int pos = 0;
  // get the velocities, accelerations, and previous velocities
  if (dsp == 0)  {
    dsp = new Vector [numSub];
    vel = new Vector [numSub];
    acc = new Vector [numSub];
    pVel = new Vector [numSub];
  }

  for(int iSub = 0; iSub < numSub; iSub++) {

    dsp[iSub].setData(tmpDisp->subData(iSub), tmpDisp->subLen(iSub));

    vel[iSub].setData(dState.getVeloc().subData(iSub), 
                      dState.getVeloc().subLen(iSub));

    acc[iSub].setData(dState.getAccel().subData(iSub), 
                      dState.getAccel().subLen(iSub));

    pVel[iSub].setData(dState.getPrevVeloc().subData(iSub), 
                       dState.getPrevVeloc().subLen(iSub));
  }

  for(int iNeigh = 0; iNeigh < nbrReceivingFromMe; iNeigh++) {

    // get fluid mpi to send disps to
    int mpiNum = idSendTo[iNeigh];

    int origPos = pos;

    // compute displacements for each match point
    for(int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;
      int locElem = sndTable[mpiNum][iData].elemNum;

      int iSub = locSub;
      State localState(cdsa[iSub], dsa[iSub], usrDefDisps[iSub],
                       usrDefVels[iSub], dsp[iSub],
                       vel[iSub], acc[iSub], pVel[iSub]);
      GeomState *geomState = (distrGeomState) ? (*distrGeomState)[iSub] : 0;

      if(!useFaceElem) {
        thisElement = (*eset[locSub])[locElem];
        thisElement->computeDisp(*cs[locSub], localState,
                                sndTable[mpiNum][iData], buffer+pos, geomState);
      }
      else {
        thisFaceElem = (*feset)[locElem];
        thisFaceElem->computeDisp(*cs[locSub], localState, sndTable[mpiNum][iData],
                                  buffer+pos, geomState, fnId2[locSub]); // TODO fnId2 must map from embedded surface to locSub node numbers
      }

      pos += 6;
    }
  
    if(tag < 0) tag = STTOFLMT+( (sndParity > 0) ? 1 : 0);

    toFluid = 1;
    for(int ip = origPos; ip < pos; ip += 6) {
      xxx += buffer[ip+0]*buffer[ip+0] + buffer[ip+1]*buffer[ip+1] +
             buffer[ip+2]*buffer[ip+2];
      yyy += buffer[ip+3]*buffer[ip+3] + buffer[ip+4]*buffer[ip+4] +
             buffer[ip+5]*buffer[ip+5];
    }

    fluidCom->sendTo(mpiNum, tag, buffer+origPos, pos-origPos);
  }
  fluidCom->waitForAllReq();
  //if(verboseFlag)
  //  filePrint(stderr, "Sending %e and %e to %d CPUs\n", xxx, yyy, nbrReceivingFromMe);
  flipSndParity();
}

//---------------------------------------------------------------------

void DistFlExchanger::sendParam(int algnum, double step, double totaltime, 
		                int rstinc, int _isCollocated, double _a[2])
{
  //int zero = 0;
  int TNd  = 0;
  int thisNode = structCom->myID();

  double buffer[5];
  buffer[0] = (double) algnum;
  buffer[1] = (algnum==5)? step/2 : step;
  //buffer[1] = step;
  buffer[2] = totaltime;
  buffer[3] = (double) rstinc;
  buffer[4] = (double) 2; // Used to be AeroScheme now always conservative
  int msglen = 5;
  int tag = 3000;

  if(thisNode == 0) {
     fluidCom->sendTo(TNd, tag, buffer, msglen);
     fluidCom->waitForAllReq();
  }

  isCollocated = _isCollocated;
  dt = step;
  if (algnum > 4) 
    alpha[0] = alpha[1] = 0;
  else  {
    alpha[0] = _a[0];
    alpha[1] = _a[1];
  }
}

//---------------------------------------------------------------------

void
DistFlExchanger::sendTempParam(int algnum, double step, double totaltime,
                               int rstinc, double alphat[2])
{
  int TNd  = 0;
  int thisNode;
  thisNode = structCom->myID();

  double buffer[5];
  buffer[0] = (double) algnum;
  buffer[1] = step;
  buffer[2] = totaltime;
  buffer[3] = (double) rstinc;
  buffer[4] = (double) 2; // Used to be AeroScheme now always conservative
  int msglen = 5;
  int tag = 3000;

// Send to SUBROUTINE GETTEMPALG of Fluid Code
  if(thisNode == 0){
   fluidCom->sendTo(TNd, tag, buffer, msglen);
   fluidCom->waitForAllReq();
   }

  dtemp = step;
  alph[0] = alphat[0];
  alph[1] = alphat[1];

//  waitOnSend();
}

//---------------------------------------------------------------------

double
DistFlExchanger::getFluidLoad(DistrVector& force, int tIndex, 
                              double time, double alphaf, int& iscollocated,
                              DistrGeomState* distrGeomState)
{
  aforce[0] = aforce[1] = aforce[2] = 0.0;

  FaceElemSet *feset;
  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
  }

  Element     *thisElement;
  FaceElement *thisFaceElem;

  for(int iNeigh = 0; iNeigh < nbrReceivingFromMe; iNeigh++) {

    int tag =  FLTOSTMT + ((rcvParity > 0) ? 1 : 0) ;
    toFluid = 1;

    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    // get fluid mpi process
    int mpiNum = rInfo.cpu;

    for(int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;  
      int locElem = sndTable[mpiNum][iData].elemNum;

      GeomState *geomState = (distrGeomState) ? (*distrGeomState)[locSub] : 0;

      int nDof;
      if(!useFaceElem) {
        thisElement = (*eset[locSub])[locElem];
        thisElement->getFlLoad(*(cs[locSub]), sndTable[mpiNum][iData], 
                               buffer+3*iData, localF, geomState);
        nDof = thisElement->numDofs();
      }
      else {
        thisFaceElem = (*feset)[locElem];
        thisFaceElem->getFlLoad(sndTable[mpiNum][iData], buffer+3*iData, localF);
        nDof = thisFaceElem->numDofs();
      }

      int *dof = sndTable[mpiNum][iData].dofs;
      for(int iDof = 0; iDof < nDof; ++iDof)
        if(dof[iDof] >= 0)
          force.subData(locSub)[dof[iDof]] += localF[iDof];
         
      aforce[0] += buffer[3*iData];
      aforce[1] += buffer[3*iData+1];
      aforce[2] += buffer[3*iData+2];
    }
  }

  flipRcvParity();

  if (oinfo) {
    if (tIndex % oinfo->interval == 0 && oinfo->filptr != NULL) {
      structCom->reduce(3,aforce);
      filePrint(oinfo->filptr,"%e   ",time);
      filePrint(oinfo->filptr,"%e %e %e\n",aforce[0],aforce[1],aforce[2]);
      fflush(oinfo->filptr);
    }
  }

  // ML & KP For 'corrected' aeroelastic force
  iscollocated = (isCollocated) ? 1 : 0;
  return (time+alphaf*dt);
}

//---------------------------------------------------------------------

void
DistFlExchanger::sendModeFreq(double *modfrq, int nummod)
{

        int TNd  = 0;

        int ThisNode = structCom->myID();

        int Type = 1200;

        double *sBuffer = new double[1+nummod];
        sBuffer[0]   = (double) nummod;

        int mod;
        for (mod = 0; mod < nummod; mod++)
                sBuffer[1+mod] = modfrq[mod];

        int Pos  = 1 + nummod;

        if (ThisNode == 0) {
                // _FORTRAN(nsedoc)(zero, Type, sBuffer, Pos,
                //                 TNd, toFluid);
           fluidCom->sendTo(TNd, Type, sBuffer, Pos);
           fluidCom->waitForAllReq();
        }

}

void
DistFlExchanger::sendModeShapes(int numModes, int numN, double (**v)[6],
                                SysState<DistrVector>& st, double ampFactor)
{
  cerr << "ERROR: DistFlExchanger::sendModeShapes is not implemented\n"; 
  exit(-1);
/*
  int iMode;
  for(iMode = 0; iMode < numModes; ++iMode) {
     // Create the state
     Vector &d = st.getDisp();
     d.zero();
     int iNode, dof;
     for(iNode = 0; iNode < numN; ++iNode) {
       dof = dsa->locate(iNode, DofSet::Xdisp);
       if(dof >= 0) d[dof] = v[iMode][iNode][0]*ampFactor;
       dof = dsa->locate(iNode, DofSet::Ydisp);
       if(dof >= 0) d[dof] = v[iMode][iNode][1]*ampFactor;
       dof = dsa->locate(iNode, DofSet::Zdisp);
       if(dof >= 0) d[dof] = v[iMode][iNode][2]*ampFactor;
       dof = dsa->locate(iNode, DofSet::Xrot);
       if(dof >= 0) d[dof] = v[iMode][iNode][3]*ampFactor;
       dof = dsa->locate(iNode, DofSet::Yrot);
       if(dof >= 0) d[dof] = v[iMode][iNode][4]*ampFactor;
       dof = dsa->locate(iNode, DofSet::Zrot);
       if(dof >= 0) d[dof] = v[iMode][iNode][5]*ampFactor;
     }
     //fprintf(stderr, "Total mode disp %e\n", d*d);
     sendDisplacements(st, 1201+iMode);
  }
*/
}

void
DistFlExchanger::thermoread(int bLen)
{
// Initialize the buffer

  buffLen = bLen;
  buff = new double[buffLen];
}


void
DistFlExchanger::sendStrucTemp(DistrVector& tempsent)
{
  // for now we assume that the structure and thermal models have the same mesh and the same decomposition
  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

// Sends temperature to mechanical structure
// for now we assume that the structure and thermal models have the same mesh and the same decomposition

   int thisNode;
   thisNode = structCom->myID();
   int tag = STTOSTMT;
   int zero = 0;

   int i;
   int pos=tempsent.size();

   for(i=0; i<pos; i++) {
     buff[i]=tempsent[i];
    }

//    if (thisNode==0)
    //fprintf(stderr,"++++ Sending TempLoad \n");

    heatStructCom->sendTo(thisNode, tag, buff, pos);

    //fprintf(stderr,"++++ Finished Sending TempLoad \n");
}

void
DistFlExchanger::getStrucTemp(double* temprcvd)
{
// Receives temperature from Thermal Elements
// temprcvd are NODAL temperatures

    int thisNode;
    thisNode = structCom->myID();
    int tag = STTOSTMT;
    int rsize;
    int fromNd;

//    if(thisNode==0){
    //fprintf(stderr,"---- getting TempLoad \n");
    RecInfo rInfo = heatStructCom->recFrom(tag, buff, buffLen);
    fromNd = rInfo.cpu;
    rsize = rInfo.len;

    //fprintf(stderr,"++++ Finished Receiving TempLoad \n");

    int i;
    for (i=0; i<rsize; i++) { temprcvd[i]=buff[i] ;
      //fprintf(stderr,"**** received  = %f %d %d \n",buff[i], buffLen, rsize);
    }

}

int
DistFlExchanger::cmdCom(int commandFlag)
{
  int returnFlag = 0;
  int FldNd  = 0;
  int tag;
  int thisNode = structCom->myID();
  double buffer[1] = { double(commandFlag) };
  int msglen = 1;

  if(thisNode == 0) {
    // fprintf(stderr,"\nSTC: sending command to Fluid: %d\n",commandFlag);   
    // fflush(stderr);
    tag = STCMDMSG;
    fluidCom->sendTo(FldNd, tag, buffer, msglen);
    tag =  FLCMDMSG;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
    returnFlag = (int) buffer[0];
    // fprintf(stderr,"\nSTC: obtained command from Fluid: %d\n",returnFlag);   
    // fflush(stderr);
  }

  return returnFlag;
}

void
DistFlExchanger::sendTemperature(SysState<DistrVector>& dState)
{
  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  if (tmpDisp == 0)
    tmpDisp = new DistrVector(dState.getDisp());

  *tmpDisp = dState.getDisp();
  tmpDisp->linAdd(dtemp*alpha[0], dState.getVeloc(), dtemp*alpha[1], dState.getPrevVeloc());

  double xxx = 0, yyy = 0;
  int pos = 0;

  // get the temperature, velocities and previous velocities
  if (dsp == 0)  {
    dsp = new Vector [numSub];
    vel = new Vector [numSub];
    pVel = new Vector [numSub];
  }

  for (int iSub = 0; iSub < numSub; iSub++)  {

    dsp[iSub].setData(tmpDisp->subData(iSub), tmpDisp->subLen(iSub));

    vel[iSub].setData(dState.getVeloc().subData(iSub), 
	      	     	   dState.getVeloc().subLen(iSub));

    pVel[iSub].setData(dState.getPrevVeloc().subData(iSub), 
		      	    dState.getPrevVeloc().subLen(iSub));
  }

  for (int iNeigh = 0; iNeigh < nbrReceivingFromMe; iNeigh++) {

    // get fluid mpi to send disps to
    int mpiNum = idSendTo[iNeigh];

    int origPos = pos;

    // compute displacements for each match point
    for (int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;
      int locElem = sndTable[mpiNum][iData].elemNum;
      Element *thisElement = (*eset[locSub])[locElem];

      int iSub = locSub;
      State localState(cdsa[iSub], dsa[iSub], (double *) 0 /*usrDefDisps[iSub]*/,
                       dsp[iSub], vel[iSub], pVel[iSub]);
      thisElement->computeTemp(*cs[locSub], localState,
      			       sndTable[mpiNum][iData].xy, buffer+pos);
      pos += 1;

    }
  
    int tag = STTOFLHEAT;

    fluidCom->sendTo(mpiNum, tag, buffer+origPos, pos-origPos);
  }
  fluidCom->waitForAllReq();
}

double
DistFlExchanger::getFluidFlux(DistrVector& flux, double time, double &bflux)
{
  for (int iNeigh = 0; iNeigh < nbrReceivingFromMe; iNeigh++) {

    int tag =  FLTOSTHEAT;
    toFluid = 1;

    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    // get fluid mpi process
    int mpiNum = rInfo.cpu;

    for (int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;
      int locElem = sndTable[mpiNum][iData].elemNum;
      Element *thisElement = (*eset[locSub])[locElem];

      thisElement->getFlFlux(sndTable[mpiNum][iData].xy, buffer+iData, localF);

      int nDof = thisElement->numDofs();
      int *dof = sndTable[mpiNum][iData].dofs;
      for (int iDof = 0; iDof < nDof; ++iDof)
        if (dof[iDof] >= 0) {
          flux.subData(locSub)[dof[iDof]] += localF[iDof];
          bflux += localF[iDof];
        }
    }
  }

  return bflux;
}

//KW: send the embedded wet surface to fluid
void
DistFlExchanger::sendEmbeddedWetSurface()
{
  if(structCom->myID()!=0)
    return; /* do nothing */
  if(!useFaceElem) {
    fprintf(stderr,"ERROR: Embedded Wet Surface undefined! Aborting...\n"); exit(-1);}

  // info about surface
  FaceElemSet& feset = surface->GetFaceElemSet();
  CoordSet&   fnodes = surface->GetNodeSet();
  int*          fnId = surface->GetPtrGlNodeIds();
  map<int,int>*  g2l = surface->GetPtrGlToLlNodeMap();
  int         nNodes = fnodes.size();
  int         nElems = surface->nFaceElements();
  bool    renumbered = surface->IsRenumbered();

  // data preparation
  int    buf[2] = {nNodes, nElems};
  double nodes[3*nNodes];
  int    elems[3*nElems];
  for(int i=0; i<nNodes; i++) {
    nodes[3*i]   = (*globalCoords)[fnId[i]]->x; // TODO do it without globalCoords
    nodes[3*i+1] = (*globalCoords)[fnId[i]]->y;
    nodes[3*i+2] = (*globalCoords)[fnId[i]]->z;
  }

  if(renumbered)
    for(int i=0; i<nElems; i++) {
      FaceElement *ele = feset[i];
      elems[3*i]   = ele->GetNode(0);
      elems[3*i+1] = ele->GetNode(1);
      elems[3*i+2] = ele->GetNode(2);
    }
  else
    for(int i=0; i<nElems; i++) {
      FaceElement *ele = feset[i];
      elems[3*i]   = (*g2l)[ele->GetNode(0)];
      elems[3*i+1] = (*g2l)[ele->GetNode(1)];
      elems[3*i+2] = (*g2l)[ele->GetNode(2)];
    }

  // send the package sizes
  fluidCom->sendTo(0, 555/*tag*/, buf, 2);
  fluidCom->waitForAllReq();

  // send the node set
  fluidCom->sendTo(0, 666/*tag*/, nodes, nNodes*3);
  fluidCom->waitForAllReq();

  // send the element set
  fluidCom->sendTo(0, 888/*tag*/, elems, nElems*3);
  fluidCom->waitForAllReq();
}

