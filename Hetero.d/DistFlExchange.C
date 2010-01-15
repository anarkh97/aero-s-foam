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

#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Driver.d/DynamProbType.h>
#include <Feti.d/DistrVector.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/SysState.h>

//--------------------------------------------------------------------

// global variables
extern Communicator *structCom, *fluidCom;

int toFluid = 1;
int fromFluid = 1;
int toStruc = 2;
int fromTemp = 3;

//---------------------------------------------------------------------

DistFlExchanger::DistFlExchanger(CoordSet **_cs, Elemset **_eset, 
	DofSetArray **_cdsa, DofSetArray **_dsa,  OutputInfo *_oinfo) {

  cs = _cs;
  cdsa   = _cdsa;
  dsa = _dsa;
  oinfo = _oinfo;
  tmpDisp = 0;

  eset = _eset;

  dsp = 0;
  vel = 0;
  acc = 0;
  pVel = 0;


}

//---------------------------------------------------------------------

MatchMap DistFlExchanger::getMatchData()  {

  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);
  int *numMatch = geoSource->getNumMatchData();

  MatchMap globMatches;

  // create global sub to local sub map
  int totSub = cpuToSub->numConnect();
  int *glToLocSubMap = new int[totSub];
 
  int iSub;
  for (iSub = 0; iSub < totSub; iSub++)
    glToLocSubMap[iSub] = -1;

  for (iSub = 0; iSub < numSub; iSub++)
    glToLocSubMap[ (*cpuToSub)[myCPU][iSub] ] = iSub;

  for (iSub = 0; iSub < numSub; iSub++)  {

    int glSub = (*cpuToSub)[myCPU][iSub];
    int locSub = glToLocSubMap[glSub];
    if (locSub < 0)
      fprintf(stderr,"*** ERROR: Bad Global to Local Sub Map\n");

    MatchData *matchData = geoSource->getMatchData(locSub);

    for (int iMatch= 0; iMatch < numMatch[locSub]; iMatch++)  {

      int gPos = matchData[iMatch].glPos;
      InterpPoint iPoint;
      iPoint.subNumber = iSub;
      iPoint.elemNum = matchData[iMatch].elemNum;
      iPoint.xy[0] = matchData[iMatch].xi;
      iPoint.xy[1] = matchData[iMatch].eta;
      globMatches[gPos] = iPoint;
    }
  }

  return globMatches;
}

//----------------------------------------------------------------------

// This routine negotiate with the fluid codes which match points go where
void DistFlExchanger::negotiate()  {

  
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

  MatchMap globMatches = getMatchData();

  // # of fluid mpi's receiving data from this mpi
  numFluidNeighbors = 0;

  int *index = new int[maxSize];
  fprintf(stderr,"struct %d allocating sndTable of size %d\n", structCom->myID(), nSender);
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
    //int sender = consOrigin[fromNd];  // gives fluid mpi #
    int sender = rInfo.cpu;
    
    //fprintf(stderr,"Struct cpu %d receives %d matches from fluid %d\n", structCom->myID(), nbSendTo[sender], sender); 

    // determine how many matches are in this mpi process
    int numHits = 0;
    int iMatch;
    for (iMatch = 0; iMatch < nbSendTo[sender]; iMatch++)
      if (globMatches.find(index[iMatch]) != globMatches.end()) 
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
        if (globMatches.find(index[iMatch]) != globMatches.end())  {
          matchList[iHit++] = iMatch;
        
        }

      // reset nbSendTo to actual number of matches in this mpi 
      nbSendTo[sender] = numHits;

      // assign match data to the sndTable
      //for (int ipt = 0; ipt < numHits; ++ipt) {
      for (int ipt = 0; ipt < nbSendTo[sender]; ++ipt) {
       
        int arrayPos = index[matchList[ipt]];
        sndTable[sender][ipt].subNumber = globMatches[arrayPos].subNumber; 
        sndTable[sender][ipt].elemNum = globMatches[arrayPos].elemNum; 
        sndTable[sender][ipt].xy[0] = globMatches[arrayPos].xy[0]; 
        sndTable[sender][ipt].xy[1]= globMatches[arrayPos].xy[1]; 

        // assign dofs to element 
        int locSub = sndTable[sender][ipt].subNumber;
        int locElem = sndTable[sender][ipt].elemNum;
        Element *thisElement = (*eset[locSub])[locElem];
        int nDof = thisElement->numDofs();
        totalNDof += nDof;
        if (nDof > maxNDof) 
          maxNDof = nDof;
      }
      numFluidNeighbors++;
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
        Element *thisElement = (*eset[locSub])[locElem];

        sndTable[sender][iData].dofs = array;
        thisElement->dofs(*(cdsa[locSub]), array);
        array += thisElement->numDofs();;

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
  idSendTo = new int[numFluidNeighbors]; 

  int count = 0;
  for (iFluid = 0; iFluid < nSender; iFluid++)
    if (sndTable[iFluid] != 0)
      idSendTo[count++] = iFluid;

  delete [] index;
  buffer = new double[6*totSize];
  setBufferLength(6*totSize);
}

//---------------------------------------------------------------------

void DistFlExchanger::sendDisplacements(SysState<DistrVector> &dState, 
		double **usrDefDisps, double **usrDefVels,  int tag)  {

  Connectivity *cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  if (tmpDisp == 0)
    tmpDisp = new DistrVector(dState.getDisp());

  *tmpDisp = dState.getDisp();
  tmpDisp->linAdd(dt*alpha[0], dState.getVeloc(), dt*alpha[1], dState.getPrevVeloc());

  //fprintf(stderr, "Disp Norm %e Veloc Norm %e\n", newState.getDisp()*newState.getDisp(), newState.getVeloc()*newState.getVeloc());

  double xxx = 0;
  int pos = 0;


  // get the velocities, accelerations, and previous velocities
  if (dsp == 0)  {
    dsp = new Vector [numSub];
    vel = new Vector [numSub];
    acc = new Vector [numSub];
    pVel = new Vector [numSub];
    newState = new State [numSub];


    newState = new State [numSub];
  }

  for (int iSub = 0; iSub < numSub; iSub++)  {

    dsp[iSub].setData(tmpDisp->subData(iSub), tmpDisp->subLen(iSub));

    vel[iSub].setData(dState.getVeloc().subData(iSub), 
	      	     	   dState.getVeloc().subLen(iSub));

    acc[iSub].setData(dState.getAccel().subData(iSub), 
		       	   dState.getAccel().subLen(iSub));

    pVel[iSub].setData(dState.getPrevVeloc().subData(iSub), 
		      	    dState.getPrevVeloc().subLen(iSub));
/*
    newState[iSub] = State(cdsa[iSub], dsa[iSub], usrDefDisps[iSub], 
	 		usrDefVels[iSub], dsp[iSub], vel[iSub], 
			acc[iSub], pVel[iSub]);
*/
  }

  for (int iNeigh = 0; iNeigh < numFluidNeighbors; iNeigh++) {

    // get fluid mpi to send disps to
    int mpiNum = idSendTo[iNeigh];

    int origPos = pos;

    // compute displacements for each match point
    for (int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;
      int locElem = sndTable[mpiNum][iData].elemNum;
      Element *thisElement = (*eset[locSub])[locElem];

      int iSub = locSub;
      State localState(cdsa[iSub], dsa[iSub], usrDefDisps[iSub],
                        usrDefVels[iSub], dsp[iSub], vel[iSub],
                        acc[iSub], pVel[iSub]);
      thisElement->computeDisp(*cs[locSub], localState,
			       sndTable[mpiNum][iData], buffer+pos);
      pos += 6;

    }
  
    if (tag < 0) 
      tag = STTOFLMT+( (sndParity > 0) ? 1 : 0);

    toFluid = 1;
    for (int ip = origPos; ip < pos; ip += 6)
      xxx += buffer[ip+3]*buffer[ip+3] + buffer[ip+4]*buffer[ip+4] +
            buffer[ip+5]*buffer[ip+5];

    fluidCom->sendTo(mpiNum, tag, buffer+origPos, pos-origPos);
  }
  fluidCom->waitForAllReq();
  fprintf(stderr, "Sending %e to %d CPUs\n",xxx, numFluidNeighbors);
  flipSndParity();

}

//---------------------------------------------------------------------

void DistFlExchanger::sendParam(int algnum, double step, double totaltime, 
		int rstinc, int _isCollocated, double _a[2])  {

  //int zero = 0;
  int TNd  = 0;
  int thisNode = structCom->myID();

  double buffer[5];
  buffer[0] = (double) algnum;
  buffer[1] = step;
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

double DistFlExchanger::getFluidLoad(DistrVector &force, int tIndex, 
			double time, double alphaf, int &iscollocated)  {

  aforce[0] = 0.0;
  aforce[1] = 0.0;
  aforce[2] = 0.0;

  //int zero = 0;
  //Connectivity *cpuToSub = geoSource->getCpuToSub();
  //int myCPU = structCom->myID();
  //int numSub = cpuToSub->num(myCPU);

  for (int iNeigh = 0; iNeigh < numFluidNeighbors; iNeigh++) {

    int tag =  FLTOSTMT + ((rcvParity > 0) ? 1 : 0) ;
    toFluid = 1;

    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    // get fluid mpi process
    int mpiNum = rInfo.cpu;

    for (int iData = 0; iData < nbSendTo[mpiNum]; ++iData) {

      int locSub = sndTable[mpiNum][iData].subNumber;  
      int locElem = sndTable[mpiNum][iData].elemNum;
      Element *thisElement = (*eset[locSub])[locElem];
      
      thisElement->getFlLoad(*(cs[locSub]), sndTable[mpiNum][iData], 
			     buffer+3*iData, localF);

      int nDof = thisElement->numDofs();
      int *dof = sndTable[mpiNum][iData].dofs;
      for (int iDof = 0; iDof < nDof; ++iDof)
        if (dof[iDof] >= 0)
          force.subData(locSub)[dof[iDof]] += localF[iDof];
         
      aforce[0] += buffer[3*iData];
      aforce[1] += buffer[3*iData+1];
      aforce[2] += buffer[3*iData+2];
    }
  }

  //fprintf(stderr,"First Buffer Norms: %e %e %e \n", aforce[0],aforce[1],aforce[2]);

  flipRcvParity();

  // KHP
  if (oinfo) {
    fprintf(oinfo->filptr,"%e   ",time);
    fprintf(oinfo->filptr,"%e %e %e\n",aforce[0],aforce[1],aforce[2]);
    fflush(oinfo->filptr);
  }


  // ML & KP For 'corrected' aeroelastic force
 iscollocated = (isCollocated) ? 1 : 0;
 return (time+alphaf*dt);
}

//---------------------------------------------------------------------

int DistFlExchanger::cmdCom(int commandFlag)
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


