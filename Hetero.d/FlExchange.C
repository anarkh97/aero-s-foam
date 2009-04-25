#include <Utils.d/dbg_alloca.h>
#include <stdlib.h> 

#include <Math.d/SparseMatrix.h>
#include <Hetero.d/FlExchange.h>
#include <Hetero.d/FilteredFile.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>

#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Comm.d/Communicator.h>

#include <Driver.d/GeoSource.h>

// double dummyVariable[10];

char *RECEIVE_LIST_KW = "RCVF";
char *SUBDOMAIN_KW = "SUBD";
char *SEND_LIST_KW = "SNDF";

extern Communicator *structCom, *fluidCom, *heatStructCom;
extern int verboseFlag;

FlExchanger::FlExchanger(Elemset& _eset, DofSetArray *_dsa, 
                         OutputInfo *_oinfo) : eset(_eset)
{
 dsa      = _dsa;
 oinfo    = _oinfo;
 tmpDisp  = 0;

}


double
FlExchanger::getFluidLoad(CoordSet &cs, Vector &force, int tIndex, double time,
                          double alphaf, int &iscollocated)
{

 aforce[0] = aforce[1] = aforce[2]=0.0;
 int i,j,iDof;

 for(i=0; i<nbrReceivingFromMe; i++) {
     int fromNd;
     int tag =  FLTOSTMT + ((rcvParity > 0) ? 1 : 0) ;
     //_FORTRAN(nrecoc)(zero, tag, buffer, bufferLen, rsize, fromNd, toFluid);
     RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
     fromNd = rInfo.cpu;
     int origin = consOrigin[fromNd];
     for(j=0; j<nbSendTo[origin]; ++j) {
        Element *thisElement = eset[sndTable[origin][j].elemNum];
        thisElement->getFlLoad(cs, sndTable[origin][j], buffer+3*j, localF);
        int nDof = thisElement->numDofs();
        int *dof = sndTable[origin][j].dofs;
        for(iDof=0; iDof<nDof; ++iDof)
          if(dof[iDof] >= 0) {
            force[dof[iDof]] += localF[iDof];
          }
        aforce[0] += buffer[3*j];
        aforce[1] += buffer[3*j+1];
        aforce[2] += buffer[3*j+2];
     }
 }
 
 flipRcvParity();


 // KHP
 if(oinfo) {
   if (tIndex % oinfo->interval == 0) {
     fprintf(oinfo->filptr,"%e   ",time);
     fprintf(oinfo->filptr,"%e %e %e\n",aforce[0],aforce[1],aforce[2]);
     fflush(oinfo->filptr);
   }
 }


 // ML & KP For 'corrected' aeroelastic force
 iscollocated = (isCollocated) ? 1 : 0;
 return (time+alphaf*dt);
}

/*
        Inform the Fluid Code About the Selected AeroeElasticity
        Algorithms and Structure Time Step
*/

extern int mflag;

void
FlExchanger::sendDisplacements( CoordSet & cs, State &state, int tag)  {

 if(tmpDisp == 0)
 tmpDisp = new Vector(state.getDisp());

 *tmpDisp = state.getDisp();
 //fprintf(stderr, "Now %e\n", *tmpDisp * *tmpDisp);
 tmpDisp->linAdd(dt*alpha[0], state.getVeloc(), dt*alpha[1], state.getPrevVeloc());
 State newState(state, *tmpDisp);
 if (verboseFlag)
   fprintf(stderr, "Disp Norm %e Veloc Norm %e\n", newState.getDisp()*newState.getDisp(), newState.getVeloc()*newState.getVeloc());

 int i,j;
 double xxx = 0;
 int pos = 0;
 for(i=0; i < nbrReceivingFromMe; i++) {

   int origPos = pos;
   //buffer[0] = mynode+1; buffer[1] = 6;
   for(j=0; j < nbSendTo[i]; ++j) {
     Element *thisElement = eset[sndTable[i][j].elemNum];

     thisElement->computeDisp(cs, newState, sndTable[i][j], buffer+pos);

     pos += 6;

   }
   
   // int tag = STTOFLMT+( (sndParity > 0) ? 1 : 0);
   if(tag < 0) tag = STTOFLMT+( (sndParity > 0) ? 1 : 0);
   int fluidNode  = idSendTo[i];

   for(int ip = origPos; ip < pos; ip += 6)
     xxx += buffer[ip+3]*buffer[ip+3] + buffer[ip+4]*buffer[ip+4] +
            buffer[ip+5]*buffer[ip+5];

   //_FORTRAN(nsedoc)(zero, tag, buffer, pos, fluidNode, toFluid);
   fluidCom->sendTo(fluidNode, tag, buffer+origPos, pos-origPos);

 }
 fluidCom->waitForAllReq();
 if (verboseFlag)
   fprintf(stderr, "Sending %e to %d CPUs\n",xxx, nbrReceivingFromMe);
 flipSndParity();

}

void
FlExchanger::sendModeShapes(int numModes, int numN, double (**v)[6],
               CoordSet &cs, State &st, double ampFactor)
{
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
     fprintf(stderr, "Total mode disp %e\n", d*d);
     sendDisplacements(cs, st, 1201+iMode);
  }
}
void
FlExchanger::waitOnSend()
{
 //_FORTRAN(nwonsd)();
 fluidCom->waitForAllReq();
}



void
FlExchanger::sendParam(int algnum, double step, double totaltime,
                       int rstinc, int _isCollocated, double _a[2])
{
  int TNd  = 0;
  int thisNode;
  //_FORTRAN(getnod)(thisNode);
  thisNode = structCom->myID();

  double buffer[5];
  buffer[0] = (double) algnum;
  buffer[1] = (algnum==5)? step/2 : step;
  //buffer[1] = step;
  buffer[2] = totaltime;
  buffer[3] = (double) rstinc;
  buffer[4] = (double) 2; // Used to be AeroScheme now always conservative
  int msglen = 5;
  int tag = 3000;

  fprintf(stderr, " ... Sending algorithm %d\n", algnum);

  if(thisNode == 0) {
     // _FORTRAN(nsedoc)(zero, tag, buffer, msglen,
     //                            TNd, toFluid);
     // _FORTRAN(nwonsd)();
     fluidCom->sendTo(TNd, tag, buffer, msglen);
     fluidCom->waitForAllReq();
  }

  isCollocated = _isCollocated;
  dt = step;
  alpha[0] = _a[0];
  alpha[1] = _a[1];
}

void
FlExchanger::sendTempParam(int algnum, double step, double totaltime,
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

  waitOnSend();
}



void FlExchanger::sendModeFreq(double *modfrq, int nummod)
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

/*
void FlExchanger::sendModeShapes(CoordinateSet & Coords,
                         FreedomSet &Frdms,
                         FixedSet &Presc,
                           StateSet &State,
                           ActiveDof &DofList,
                           int nummod,
                           int numnod,
                           double *ModShp)
{
 int zero = 0;
 int mynode;
 _FORTRAN(getnod)(mynode);
 double veloc[1000];

  int ndf = 6;
        for (unsigned int mod = 0; mod < nummod; mod++)
                {

                double* CurrModShp = ModShp + numnod * ndf * mod;
                State.FillModeShape(CurrModShp,
                                    Frdms.numdof(),
                                    Frdms.idm());


 int i,j;
 for(i=0; i < nbrReceivingFromMe; i++) {
   buffer[0] = mynode+1; buffer[1] = 6;
   int pos = 2;
   for(j=0; j < nbSendTo[i]; ++j) {
     double vp[3];
     Element *thisElement = allElements[sndTable[i][j].elemNum];
     thisElement->Locate(Coords,Frdms,Presc,State, DofList);
     thisElement->SetVeloc(Frdms,State, DofList,veloc);
     thisElement->ComputeDisp(  veloc,
                                sndTable[i][j].xy,
                                buffer+pos, vp);
     buffer[pos] += dt*(alpha[0]*buffer[pos+3]+ alpha[1]*vp[0]);
     buffer[pos+1] += dt*(alpha[0]*buffer[pos+4]+ alpha[1]*vp[1]);
     buffer[pos+2] += dt*(alpha[0]*buffer[pos+5]+ alpha[1]*vp[2]);
     pos += 6;
   }
   int tag = 1201+mod;
   int fluidNode  = idSendTo[i];
   // fprintf(stderr," STRUCT : Ready to send to Node %4d Tag %8d\n",
   //         fluidNode,tag);

  // XML TO DO
   _FORTRAN(nsedoc)(zero, tag, buffer, pos, fluidNode, toFluid);

   // fprintf(stderr," STRUCT : Done of sending to Node %4d Tag %8d\n",
   //         fluidNode,tag);
 }
}
  // XML TO DO
 _FORTRAN(nwonsd)();
}
*/

void
FlExchanger::read(int myNode, char* inputFileName)
{

 FilteredFile ffile(inputFileName);
 char *cLine;

 int j, sender, receiver, one=1;

 // Locate where the data pertaining to this node starts
 while(1) {
   if(ffile.findToken(SUBDOMAIN_KW) < 0) {
       fprintf(stderr,"Token %s not found. Aborting\n",SUBDOMAIN_KW);
       exit(one);
   }
   cLine = ffile.getLineAfterToken();
   fflush(stdout);
   int nodeNum;
   sscanf(cLine,"%d",&nodeNum);
   if(nodeNum - 1 == myNode) break;
   ffile.findToken("END");
 } 

 // Let's look for the receive list
 if(ffile.findToken(RECEIVE_LIST_KW) < 0) {
     fprintf(stderr,"Token %s not found. Aborting\n",RECEIVE_LIST_KW);
     exit(one);
 }
 cLine = ffile.getLineAfterToken();
 int numSnd;
 sscanf(cLine,"%d",&numSnd);

 int *rcvcomid  = new int[numSnd];
 int *rcvcomlen = new int[numSnd];
 int ** rcvEleList = new int*[numSnd];
 int ** rcvGaussList = new int*[numSnd];
 int actualSenders = 0; // Actual number of fluid nodes sending to me
 int maxSender = 0; // number associated with the highest actual sender
 int maxPRec = 0;   // number associated with the highest fluid process
     //   That receives from us
 int maxrec=0;      // number of the highest element receiving a pressure
 for(sender=0; sender < numSnd; ++sender)
  {
   cLine = ffile.getLine();
   sscanf(cLine, "%d%d",rcvcomid+sender,rcvcomlen+sender);

   if(rcvcomlen[sender] > 0)
    {
     rcvEleList[sender]   = new int[rcvcomlen[sender]];
     rcvGaussList[sender] = new int[rcvcomlen[sender]];
     ++actualSenders;
     if(rcvcomid[sender] > maxSender) maxSender = rcvcomid[sender];
    }
   else {
     rcvEleList[sender]   = 0;
     rcvGaussList[sender] = 0;
   }
   for(j=0; j < rcvcomlen[sender]; ++j) {
     cLine = ffile.getLine();
     sscanf(cLine,"%d%d",rcvEleList[sender]+j, rcvGaussList[sender]+j);
     if (rcvEleList[sender][j] > maxrec) maxrec = rcvEleList[sender][j];
   }
  }

 // Let's look for the send list
 if(ffile.findToken(SEND_LIST_KW) < 0) {
     fprintf(stderr,"Token %s not found. Aborting\n",SEND_LIST_KW);
     exit(one);
 }
 cLine = ffile.getLineAfterToken();
 int numRcv;
 sscanf(cLine,"%d",&numRcv);
 int *sndcomid  = new int[numRcv];
 int *sndcomlen = new int[numRcv];
 int actualReceivers = 0;
 InterpPoint **interpPoints  = new InterpPoint *[numRcv];

 // fprintf(stderr, " *** WARNING: the structure does not read the gap\n");

 for(receiver = 0; receiver < numRcv; ++ receiver) {
     cLine = ffile.getLine();
     sscanf(cLine, "%d%d",sndcomid+receiver, sndcomlen+receiver);
     if(sndcomid[receiver] > maxPRec) maxPRec = sndcomid[receiver];
     if(sndcomlen[receiver] > 0) {
       ++actualReceivers;
       interpPoints[receiver] = new InterpPoint[sndcomlen[receiver]];
     }
     else {
       interpPoints[receiver] = 0;
     }
     for(j=0; j<sndcomlen[receiver]; ++j) {
        cLine = ffile.getLine();
         //sscanf(cLine,"%d%lf%lf", &interpPoints[receiver][j].elemNum,
           //interpPoints[receiver][j].xy, interpPoints[receiver][j].xy+1);
         sscanf(cLine,"%d%lf%lf%lf%lf%lf", &interpPoints[receiver][j].elemNum,
          interpPoints[receiver][j].xy, interpPoints[receiver][j].xy+1,
               interpPoints[receiver][j].gap, interpPoints[receiver][j].gap+1,
               interpPoints[receiver][j].gap+2);
         interpPoints[receiver][j].elemNum -= 1; // Numbering convetion starts
            // at 0.

         // map global elem number to packed elemset numbering 
         interpPoints[receiver][j].elemNum = geoSource->glToPackElem(interpPoints[receiver][j].elemNum);
     }
 }


// Now create the compacted tables

//  ***  First the receive table is compacted to the actual receivers
//  *** and we create the element wet mask
 nbrSendingToMe = actualSenders;
 senderId = new int[maxSender];
 numWetElements = 0;

 bufferLen = 0;
 
 int *wetMask = new int[maxrec];
 for(j=0; j < maxrec; ++j) wetMask[j] = 0;
 int *table = new int[maxrec];
 nbGaussPoints = new int[maxrec];

 GP_Table = new int*[actualSenders];
 nbData = new int[actualSenders];

 int sId = 0;
 for(sender=0; sender < numSnd; ++sender) {

     if(rcvcomlen[sender] > 0) {
        nbData[sId] = rcvcomlen[sender];
        if(bufferLen < 2+NBPRESSDATAMAX*nbData[sId])
           bufferLen = 2+NBPRESSDATAMAX*nbData[sId];
        GP_Table[sId] = new int[rcvcomlen[sender]];
        senderId[rcvcomid[sender]-1] = sId;
        for(j = 0; j < rcvcomlen[sender]; ++j) {
           int thisEleId = rcvEleList[sender][j]-1;
           if(wetMask[thisEleId] == 0) { //We found a new wet element
               wetMask[thisEleId] = 1;
               table[thisEleId] = numWetElements;
               nbGaussPoints[numWetElements] = 0;
               ++numWetElements; 
           }
           nbGaussPoints[table[thisEleId]] ++;
        }
        sId++;
     }
  }


 pressureIndexForElem = new int[numWetElements+1];
// Now fill the tables
 pressureIndexForElem[0] = 0;
 for(j = 0; j < numWetElements; ++j)
   pressureIndexForElem[j+1] = pressureIndexForElem[j]+nbGaussPoints[j];

 sId = 0;
 for(sender=0; sender < numSnd; ++sender) {

     if(rcvcomlen[sender] > 0) {
        for(j = 0; j < rcvcomlen[sender]; ++j) {
           int globEleId = rcvEleList[sender][j]-1;
           int wetEleId = table[globEleId];
           if(rcvGaussList[sender][j] > nbGaussPoints[wetEleId])
               rcvGaussList[sender][j] = nbGaussPoints[wetEleId];
           GP_Table[sId][j] = pressureIndexForElem[wetEleId]+
                  rcvGaussList[sender][j]-1;
        }
        sId++;
     } 
  }

// *** Make the send tables
  
 nbrReceivingFromMe = actualReceivers;
 sndTable = new InterpPoint*[nbrReceivingFromMe];
 idSendTo = new int[nbrReceivingFromMe];
 nbSendTo = new int[nbrReceivingFromMe];
 // fprintf(stderr, " SIZE OF nbSendTo is %d\n",nbrReceivingFromMe);
 int rId = 0;
 int mysize = 0;
 consOrigin = new int[maxPRec+1];
 for(receiver =0; receiver < numRcv; ++receiver) {
   mysize += sndcomlen[receiver];
   if(sndcomlen[receiver] > 0) {
     if(bufferLen < 2+6*sndcomlen[receiver])
        bufferLen = 2+6*sndcomlen[receiver];
     sndTable[rId] = interpPoints[receiver];
     nbSendTo[rId] = sndcomlen[receiver];
     // fprintf(stderr," *** nbSendTo %8d = %8d\n",rId, sndcomlen[receiver]);
     idSendTo[rId] = sndcomid[receiver]-1;
     consOrigin[sndcomid[receiver]-1] = rId;
     rId++;
   }
 }

 //PHG buffer = new double[bufferLen];
 buffer = new double[mysize*6];
 pArray = new double[pressureIndexForElem[numWetElements]*NBPRESSDATAMAX];
 for(j = 0; j < pressureIndexForElem[numWetElements]*NBPRESSDATAMAX; ++j)
   pArray[j] = 0;



// Last step: Create the dof nDof arrays and allocate localF

  int totalNDof = 0;
  int maxNDof = 0;

  int i;
  for(i =0; i < nbrReceivingFromMe; i++) {
    for(j=0; j < nbSendTo[i]; ++j) {
      Element *thisElement = eset[sndTable[i][j].elemNum];
      int nDof = thisElement->numDofs();
      totalNDof += nDof;
      if(nDof > maxNDof) maxNDof = nDof;
    }
  }
  localF = new double[maxNDof];
  int *array = new int[totalNDof];
  for(i =0; i < nbrReceivingFromMe; i++) {
    for(j = 0; j < nbSendTo[i]; ++j) {
      Element *thisElement = eset[sndTable[i][j].elemNum];
      int nDof = thisElement->numDofs();
      sndTable[i][j].dofs = array;
      thisElement->dofs(*dsa,array);
      array += nDof;
    }
  }
}


// Temperature routines

void
FlExchanger::sendTemperature( CoordSet & cs,
                                State &state )
{
/*  nbrReceivingFromMe = number of fluid subdomains to which a structural
                         subdomain has to send values to.
    nbSendTo[i] = number of fluid nodes of fluid subdomain i communicating
                  to a structural subdomain
    buffer[0]   = fluid subdomain number
    buffer[1]   = number of information contained in one node  */

 //int mynode = 0;

 // Prediction

 if(tmpDisp == 0)
 tmpDisp = new Vector(state.getDisp());

 *tmpDisp = state.getDisp();
 tmpDisp->linAdd(dtemp*alph[0], state.getVeloc(), dtemp*alph[1], state.getPrevVeloc());
 State newState(state, *tmpDisp);

 int i,j;

 for(i=0; i < nbrReceivingFromMe; i++) {

   int pos = 0;

   for(j=0; j < nbSendTo[i]; ++j) {
//    fprintf(stderr,"I = %d and nbSendTo[I] = %d\n", j, nbSendTo[i]);
     Element *thisElement = eset[sndTable[i][j].elemNum];
     thisElement->computeTemp(cs, newState, sndTable[i][j].xy, buffer+pos);

//     fprintf(stderr, "Temperature Sent : %14.5e\n", buffer[pos]);
//     fprintf(stderr, "DTEMP2 : %14.5e\n", buffer[pos+1]);
     pos += 1;
   }

   int tag = STTOFLHEAT;
   int fluidNode  = idSendTo[i];

// fprintf(stderr," STRUCT : Ready to send to Node %4d\n", tag);
// for (int xyz=0; xyz < pos; xyz++)
//      printf("Sending %4d = %14.7e\n",xyz,buffer[xyz]);

   //_FORTRAN(nsedoc)(zero, tag, buffer, pos, fluidNode, toFluid);
  fluidCom->sendTo(fluidNode, tag, buffer, pos);

//  fprintf(stderr," STRUCT : Done sending to Node %4d, Buffer %d\n", tag,pos);
//  fflush(stderr);

 }
 waitOnSend();
}

void
FlExchanger::sendStrucTemp(Vector &tempsent)
{
// Sends temperature to mechanical structure
// BE CAREFUL : a reference (&) changes the value everywhere!!!

   int thisNode;
   thisNode = structCom->myID();
   int tag = STTOSTMT;
   int zero = 0;

   int i;
   int pos=tempsent.size();

   for(i=0; i<pos; i++) {
     buff[i]=tempsent[i];
    }

//    if (thisNode==0){
    //fprintf(stderr,"++++ Sending TempLoad \n");

    heatStructCom->sendTo(zero, tag, buff, pos);

    //fprintf(stderr,"++++ Finished Sending TempLoad \n");
}

void
FlExchanger::getStrucTemp(double *temprcvd)
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



double
FlExchanger::getFluidFlux(Vector &flux, double time, double &bflux)
{
/* sndTable's extensions dofs, xy, elemNum are declared in FlExchange.h
 in structure InterpPoint.
 nbrReceivingFromMe =  number of fluid subdomains */


 aflux = 0.;    /*assembled flux */
 bflux = 0.;
 int i, j, iDof;

 for(i=0; i<nbrReceivingFromMe; i++) {
     int fromNd;
     int tag =  FLTOSTHEAT ;
     int rsize;
     RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
     fromNd = rInfo.cpu;
     rsize = rInfo.len;

     int origin = consOrigin[fromNd];

// Loop Over wet points of each fluid subdomain

     for(j=0; j<nbSendTo[origin]; ++j) {
        Element *thisElement = eset[sndTable[origin][j].elemNum];
        thisElement->getFlFlux(sndTable[origin][j].xy, buffer+j, localF);
        int nDof = thisElement->numDofs();
        int *dof = sndTable[origin][j].dofs;

// Summing the fluxes in each structural wet node 

        for(iDof=0; iDof<nDof; ++iDof)
          if(dof[iDof] >= 0) {
            flux[dof[iDof]] += localF[iDof];
            bflux += localF[iDof];
          }

/* Summing the fluid fluxes. Because of conservation, the sum below should
 be equal to the sum of the assembled fluxes, ie, equal to bflux=Sum(local[iDof])
 above. Check it. */
//        fprintf(stderr, "Summing flux %d\n",j+1);
//        aflux += buffer[j];
     }
// fprintf(stderr," STRUCT : done receiving fluxes %4d, Buffer %d\n", tag,rsize);

 }

/* fprintf(stderr, "SUM of fluxes on struct %f\n ", bflux);
 fprintf(stderr, " SUM of FlDomain %d %f\n: ", i, aflux); */


/* KHP: FlExchange is called from Dynam.C. The lines above are
   for debugging purposes only */
/* if(oinfo) {
   fprintf(oinfo->filptr,"%e   ",time);
   fprintf(oinfo->filptr,"%e\n",aflux);
   fprintf(oinfo->filptr,"%e\n",bflux);
   fflush(oinfo->filptr);
 } */

 return bflux;
}


void
FlExchanger::thermoread(int &bLen)
{
// Initialize the buffer

  buffLen = bLen;
  buff = new double[buffLen];
}

void
FlExchanger::printreceiving()
{
 fprintf(stderr," ::::::: nbrReceivingFromMe = %d\n",nbrReceivingFromMe);
}

int
FlExchanger::cmdCom( int commandFlag )
{
  int returnFlag = 0;
  int FldNd  = 0;
 
  //int rsize;
  int tag;
 
  int thisNode;
  //_FORTRAN(getnod)(thisNode);

 thisNode = structCom->myID();

  double buffer[1];
  buffer[0] = (double) commandFlag;
  int msglen = 1;
  
  if(thisNode == 0) {
 
    // fprintf(stderr,"\nSTC: sending command to Fluid: %d\n",commandFlag);   
    // fflush(stderr);
     tag = STCMDMSG;  
   //_FORTRAN(nsedoc)(zero, tag, buffer,msglen,FldNd, toFluid);
   fluidCom->sendTo(FldNd, tag, buffer, msglen);
    tag =  FLCMDMSG;
   //_FORTRAN(nrecoc)(zero, tag, buffer,msglen,rsize, FldNd, toFluid);
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
    returnFlag = (int) buffer[0];
    // fprintf(stderr,"\nSTC: obtained command from Fluid: %d\n",returnFlag);   
    // fflush(stderr);
  }
  
  return returnFlag;
}

int
FlExchanger::cmdComHeat( int commandFlag )
{
  int returnFlag = 0;
  int FldNd  = 0;

  //int rsize;
  int thisNode;
  int tag;

  thisNode = structCom->myID();

  double buffer[1];
  buffer[0] = (double) commandFlag;
  int msglen = 1;
 
  if(thisNode == 0) {

   //fprintf(stderr,"\nSTC: sending command to Fluid: %d\n",commandFlag);

   tag = STCMDMSG;
   fluidCom->sendTo(FldNd, tag, buffer, msglen);
   fluidCom->waitForAllReq();

    //fprintf(stderr,"\nSTC: command to Fluid sent: %d\n",commandFlag);
    //fflush(stderr);

   tag =  FLCMDMSG;
   RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
   returnFlag = (int) buffer[0];
    // fprintf(stderr,"\nSTC: obtained command from Fluid: %d\n",returnFlag);
    // fflush(stderr);
  }

  return returnFlag;
}




// This routine negotiate with the fluid codes which match points go where
void
FlExchanger::negotiate()
{

  int totSize = 0;
  int numFl = 0;
  // _FORTRAN(hetsize)(toFluid, numFl);
  numFl = fluidCom->remoteSize();
  int iFluid;
  int *flSize = new int[numFl];
  int totmatch = 0;
  for(iFluid = 0; iFluid < numFl; ++iFluid) {
  fflush(stderr);
    int tag = FL_NEGOT;
    // double nFlMatched;
    int nFlMatched;
    int bufferLen = 1;
    int /*rsize,*/ fromNd;
    RecInfo rInfo = fluidCom->recFrom(tag, &nFlMatched, bufferLen);
    fromNd = rInfo.cpu;
    flSize[fromNd] = nFlMatched;
    //flSize[fromNd] = int(nFlMatched);
    totmatch += nFlMatched;
  }

  //PHG
  if (totmatch == 0) {
    fprintf(stderr, " *** WARNING: by-passing negotiate step\n");
    fflush(stderr);
    return;
  }

  consOrigin = new int[numFl];
  nbSendTo = new int[numFl]; // should be actual number of senders
  idSendTo = new int[numFl]; // should be actual number of senders
  int nSender = 0;
  for(iFluid = 0; iFluid < numFl; ++iFluid) {
    if(flSize[iFluid] > 0) {
      idSendTo[nSender] = iFluid;
      nbSendTo[nSender] = flSize[iFluid];
      consOrigin[iFluid] = nSender;
      totSize += flSize[iFluid];
      nSender++;
    }
  }
  nbrSendingToMe = nbrReceivingFromMe = nSender;
  //double *index = new double[totSize];
  int *index = new int[totSize];
  InterpPoint **allPt = sndTable;

  sndTable = new InterpPoint *[nSender];
  for(iFluid = 0; iFluid < nSender; ++iFluid) {
    int tag = FL_NEGOT+1;
    int bufferLen = totSize;
    int rsize, fromNd;
    //_FORTRAN(nrecoc)(zero, tag, index, bufferLen, rsize, fromNd, toFluid);
    RecInfo rInfo = fluidCom->recFrom(tag, index, bufferLen);
    fromNd = rInfo.cpu;
    rsize = rInfo.len;
    int sender = consOrigin[fromNd];
    sndTable[sender] = new InterpPoint[nbSendTo[sender]];
    for(int ipt = 0; ipt < nbSendTo[sender]; ++ipt) {
      sndTable[sender][ipt] = allPt[0][index[ipt]];
      index[ipt] = ipt;
      //sndTable[sender][ipt] = allPt[0][int(index[ipt])];
      //index[ipt] = double(ipt);
    }
    tag = FL_NEGOT;
    //double num = nbSendTo[sender];
    int num = nbSendTo[sender];
    int one = 1;
    //_FORTRAN(nsedoc)(zero, tag, &num, one, fromNd, toFluid);
    fluidCom->sendTo(fromNd, tag, &num, one);
    tag = FL_NEGOT+1;
    //_FORTRAN(nsedoc)(zero, tag, index, nbSendTo[sender], fromNd, toFluid);
    fluidCom->sendTo(fromNd, tag, index, nbSendTo[sender]);
    // To make sure we can reuse the buffers
    fluidCom->waitForAllReq();
  }
 delete [] index;
 delete [] buffer;
 buffer = new double[6*totSize];
}
