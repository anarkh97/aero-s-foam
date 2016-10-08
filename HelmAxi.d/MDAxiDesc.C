#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>

#include <Utils.d/Connectivity.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/PolygonSet.h>
#include <Element.d/Element.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Threads.d/Paral.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/MatrixTimers.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <HelmAxi.d/AxiHElem.h>
#include <HelmAxi.d/coefFourier.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <HelmAxi.d/Line2AxiSommer.h>
#include <HelmAxi.d/ScatterData.h>
#include <HelmAxi.d/MDAxiData.h>
#include <HelmAxi.d/MPCData.h>
#include <HelmAxi.d/MDAxiDesc.h>
#include <HelmAxi.d/FetiHAxi.d/FetiHAxi.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>


#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
extern ulock_t allocLock;
#endif


extern const char* problemTypeMessage[];


MDAxiDesc::MDAxiDesc(Domain *d, FourierHelmBCs *fbcs, MPCData *mpcs,
           ScatterData *scatter) {

 domain = d; 

 glBCs  = fbcs;
 glScatter = scatter;

 times  = new StaticTimers;

 glMPCs = mpcs;
 mpcToNode = 0;
 mpcToSub  = 0;

}


void 
MDAxiDesc::preProcess() {

 fprintf(stderr, "%s",problemTypeMessage[domain->probType()]);
 fprintf(stderr, " ... FETIH Solver is Selected       ...\n");

 MatrixTimers &mt = domain->getTimers();

 fprintf(stderr, " ... Reading decomposition          ...\n");
 startTimerMemory(mt.readDecomp, mt.memoryDecomp);

 // Get Decomposition File pointer
 FILE *f = geoSource->getCheckFileInfo()->decPtr;
 // If decomposition file is not specified, then
 // open default file name: DECOMPOSITION
 if (f == 0)
   f = fopen("DECOMPOSITION","r");
 getDecomp(f);

 stopTimerMemory(mt.readDecomp, mt.memoryDecomp);
 fprintf(stderr, " ... Done reading decomposition     ...\n");

 fprintf(stderr, " ... Making Subdomains              ...\n");
 makeSubD();
 getDomainConnect();
 fprintf(stderr, " ... Made Subdomains                ...\n");

 fprintf(stderr, " ... Distributing Complex BCs       ...\n");
 startTimerMemory(mt.distributeBCs, mt.memoryDistBC);
 execParal(numSub, this, &MDAxiDesc::distributeBCHelm);
 stopTimerMemory(mt.distributeBCs, mt.memoryDistBC);

 // Distribute MPC MUST be done after distributeBC
 if (glMPCs!=0) {
   startTimerMemory(mt.distributeMPCs, mt.memoryDistMPC);
   distributeMPCs();
   stopTimerMemory(mt.distributeMPCs, mt.memoryDistMPC);
 }

 //RENUMBERING
 paralApplyToAll(numSub, localData, &MDAxiData::renumberData);

 fprintf(stderr," ... Making Interface               ...\n");
 makeInterface();

 // After preparing the coarse data in each subdomain, 
 // if need be, introduction of the MPC data
 if (glMPCs!=0) {
   int nP = (numSub>glBCs->numModes) ? numSub : 2*glBCs->numModes+1;
   times->memoryRhs -= memoryUsed(); 
   paralApplyToAll(numSub, localData, &MDAxiData::makeCKCt);
   times->memoryRhs += memoryUsed(); 
   startTimerMemory(mt.createDofs, mt.memoryForm);
   execParal(nP, this, &MDAxiDesc::prepareMPCSet);
   stopTimerMemory(mt.createDofs, mt.memoryForm);
 }

 fprintf(stderr," ... Done Interface                 ...\n");
 
 fprintf(stderr, " ... Building Interface Polygons    ...\n");
 buildInterfacePolygons();
 fprintf(stderr, " ... Done Interface Polygons        ...\n");

 fprintf(stderr," ... Getting Interface Signs        ...\n");
 getInterfSigns();
 fprintf(stderr," ... Done Interface Signs           ...\n");

}


void MDAxiDesc::getDecomp(FILE *f) {

 int NbEle = domain->numElements();

 // Allocate memory for subdomain to element connectivity
 int *connect = new int[NbEle];

 // Get the number of subdomains in Decomposition file
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

 for(isub=0; isub < numSub; ++isub) {
   int nele;
   int toto = fscanf(f,"%d",&nele);
   cx[isub] = curEle;
   if(curEle + nele > NbEle) {
      fprintf(stderr," *** ERROR: This decomposition contains more "
                     "elements than the original mesh\n");
      exit(1);
   }
   int iele;
   for(iele=0; iele < nele; ++iele) {
     int toto = fscanf(f,"%d",connect+curEle);
     connect[curEle] -= 1;
     curEle++;
   }
 }

 cx[numSub] = curEle;
 subToElem = new Connectivity(numSub,cx,connect);

}


void
MDAxiDesc::makeSubD() {

 MatrixTimers &mt = domain->getTimers();
 startTimerMemory(mt.makeSubDomains, mt.memorySubdomain);
 
 mt.memoryElemToNode -= memoryUsed(); 
 Connectivity *tc = new Connectivity(domain->getEset());
 mt.memoryElemToNode += memoryUsed();

 mt.memorySubToNode -= memoryUsed();
 subToNode = subToElem->transcon(tc);
 mt.memorySubToNode += memoryUsed();

 mt.memoryNodeToSub -= memoryUsed();
 nodeToSub = subToNode->reverse();
 mt.memoryNodeToSub += memoryUsed();

 subToSub  = subToNode->transcon(nodeToSub);

 if (glMPCs) {
    mt.memoryMPCToNode -= memoryUsed();
    buildMpcToNode();
    mt.memoryMPCToNode += memoryUsed();

    mt.memoryMPCToSub -= memoryUsed();
    mpcToSub = mpcToNode->transcon(nodeToSub);
    mt.memoryMPCToSub += memoryUsed();
 }

 domain->setNumnodes( nodeToSub->csize() ); 

 localData = new MDAxiData *[numSub];

 execParal(numSub, this, &MDAxiDesc::constructSubDomains,
                   localData, domain, subToElem, subToNode);

 delete tc;
 tc = 0;

 stopTimerMemory(mt.makeSubDomains, mt.memorySubdomain);

}


void
MDAxiDesc::buildMpcToNode() {

 int i, j;
 glMPCs->totalNumMPC = 0;

 // First compute the total number of MPCs
 for (i=0; i<glMPCs->numMPCSet; ++i) {
    if (glMPCs->term[i].angle1 == glMPCs->term[i].angle2) 
        glMPCs->term[i].slices = 0; 
    glMPCs->totalNumMPC += glMPCs->term[i].slices + 1;
 }

 int Fourier;
 int totalFourier = 2*glBCs->numModes+1;        

 int *pointer = new int[glMPCs->totalNumMPC+1];
 int *target  = new int[glMPCs->totalNumMPC*totalFourier]; 
 int countMPC = 0;
 int count = 0;

 // Number the MPCs globally
 // Each MPC for a same glMPCs->term[.] are put one after the other 
 for (i=0; i<glMPCs->numMPCSet; ++i) {
    for (j=0; j<glMPCs->term[i].slices + 1; ++j) {
       pointer[countMPC++] = count;
       for (Fourier=0; Fourier<totalFourier; ++Fourier) {
           target[count++] = glMPCs->term[i].nnum;
       } 
    } 
 }

 pointer[glMPCs->totalNumMPC] = glMPCs->totalNumMPC*totalFourier;

 mpcToNode = new Connectivity(glMPCs->totalNumMPC,pointer,target);

}


void
MDAxiDesc::constructSubDomains(int iSub, MDAxiData **ld, Domain *d,
                          Connectivity *cn, Connectivity *sToN) {

 ld[iSub] = new MDAxiData(d, iSub, cn, sToN);

}


void
MDAxiDesc::getDomainConnect() {

  MatrixTimers &mt = domain->getTimers();
  startTimerMemory(mt.makeConnectivity, mt.memoryConnect);

  // Create each subdomain's interface lists
  int iSub, jSub, subJ, iNode;
  int totConnect;
  int *nConnect = (int *) dbg_alloca(sizeof(int)*numSub);
  int *flag = (int *) dbg_alloca(sizeof(int)*numSub);
  for(iSub=0; iSub < numSub; ++iSub) {
     flag[iSub] = -1;
     nConnect[iSub] = 0;
  }

  // Count connectivity
  totConnect = 0;
  for(iSub=0; iSub < numSub; ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      // Loop over the nodes connected to the current subdomain
      int thisNode = (*subToNode)[iSub][iNode];
      for(jSub = 0; jSub < nodeToSub->num(thisNode); ++jSub)
         // Loop over the subdomains connected to this Node
         if( (subJ=(*nodeToSub)[thisNode][jSub]) > iSub)
              // We only deal with connection to highered numbered 
              // subdomains so that we guarantee symmetry of lists
              {
                 if(flag[subJ] != iSub) {
                   flag[subJ] = iSub;
                   nConnect[iSub]++;
                   nConnect[subJ]++;
                   totConnect += 2;
                 }
              }
    }
  }

  int **nodeCount = new int*[numSub];
  int **connectedDomain = new int*[numSub];
  int **remoteID = new int*[numSub];

  // Allocate memory for list of connected subdomains
  int *allAlloc = new int[3*totConnect];
  for(iSub=0; iSub < numSub; ++iSub) {
     connectedDomain[iSub] = allAlloc;
     allAlloc += nConnect[iSub];
     remoteID[iSub] = allAlloc;
     allAlloc += nConnect[iSub];
     nodeCount[iSub] = allAlloc;
     allAlloc += nConnect[iSub];
  }
    
  for(iSub=0; iSub < numSub; ++iSub) {
     flag[iSub] = -1;
     nodeCount[iSub] = new int[nConnect[iSub]];
     nConnect[iSub] = 0;
  }

  int *whichLocal  = (int *) dbg_alloca(sizeof(int)*numSub);
  int *whichRemote = (int *) dbg_alloca(sizeof(int)*numSub);

  for(iSub=0; iSub < numSub; ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      int nd = (*subToNode)[iSub][iNode];
      for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub)
         if( (subJ=(*nodeToSub)[nd][jSub]) > iSub) { 
             // We only deal with connection to higher numbered subs
                if(flag[subJ] != iSub) { 
                   // Attribute location for this sub
                   flag[subJ] = iSub;
                   connectedDomain[subJ][nConnect[subJ]] = iSub;
                   connectedDomain[iSub][nConnect[iSub]] = subJ;
                   remoteID[subJ][nConnect[subJ]] = nConnect[iSub];
                   remoteID[iSub][nConnect[iSub]] = nConnect[subJ];
                   whichLocal[subJ] = nConnect[iSub]++;
                   whichRemote[subJ] = nConnect[subJ]++;
                   nodeCount[iSub][whichLocal[subJ]]=1;
                   nodeCount[subJ][whichRemote[subJ]]=1;
                } else {
                   nodeCount[iSub][whichLocal[subJ]]++;
                   nodeCount[subJ][whichRemote[subJ]]++;
                }
         }
    }
  }

  // allocate memory for interface node lists
  Connectivity **interfNode;
  interfNode = new Connectivity *[numSub];
  for(iSub=0; iSub < numSub; ++iSub)
    interfNode[iSub] = new Connectivity(nConnect[iSub], nodeCount[iSub]);

  // fill the lists
  for(iSub=0; iSub < numSub; ++iSub) {
      flag[iSub]     = -1;
      nConnect[iSub] =  0;
  }

  for(iSub=0; iSub < numSub; ++iSub) {
     for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
       int nd = (*subToNode)[iSub][iNode];
       for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub)
          if( (subJ=(*nodeToSub)[nd][jSub]) > iSub) { 
             // We only deal with connection to highered numbered subs
                if(flag[subJ] != iSub) { 
                   // Attribute location for this sub
                   flag[subJ] = iSub;
                   whichLocal[subJ] = nConnect[iSub]++;
                   whichRemote[subJ] = nConnect[subJ]++;
                   (*interfNode[iSub])[whichLocal[subJ]][0] = nd;
                   (*interfNode[subJ])[whichRemote[subJ]][0] = nd;
                   nodeCount[iSub][whichLocal[subJ]]=1;
                   nodeCount[subJ][whichRemote[subJ]]=1;
                } else {
                   int il = nodeCount[iSub][whichLocal[subJ]]++;
                   int jl = nodeCount[subJ][whichRemote[subJ]]++;
                   (*interfNode[iSub])[whichLocal[subJ]][il]  = nd;
                   (*interfNode[subJ])[whichRemote[subJ]][jl] = nd;
                }
          }
     }
  }

  for(iSub=0; iSub < numSub; ++iSub) {
    SubDomain **subds = new SubDomain *[nConnect[iSub]];

    for(jSub=0; jSub < nConnect[iSub]; ++jSub)
       subds[jSub] = localData[ connectedDomain[iSub][jSub] ];

    SComm *sc = new SComm(nConnect[iSub], connectedDomain[iSub], subds,
                         remoteID[iSub], interfNode[iSub]);

    localData[iSub]->setSComm(sc);
  }

  stopTimerMemory(mt.makeConnectivity, mt.memoryConnect);

}


void
MDAxiDesc::distributeBCHelm(int iSub) {

 localData[iSub]->extractBCs(glBCs, nodeToSub);

}


void
MDAxiDesc::distributeMPCs() {

 fprintf(stderr, " ... Distributing MPCs              ...\n");
 execParal(numSub, this, &MDAxiDesc::extractSubDomainMPCs);

}


void
MDAxiDesc::extractSubDomainMPCs(int iSub) {

 localData[iSub]->extractMPCs(glMPCs, mpcToNode, nodeToSub);

}


void 
MDAxiDesc::buildInterfacePolygons() {

 MatrixTimers &mt = domain->getTimers();
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 PolygonSet ***allPs = new PolygonSet**[numSub];

 paralApplyToAll(numSub, localData, &MDAxiData::getInterface, allPs);
 paralApplyToAll(numSub, localData, &MDAxiData::finishInterface, allPs);

 delete[] allPs;
 allPs = 0;

 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}


void
MDAxiDesc::makeInterface() {

 MatrixTimers &mt = domain->getTimers();
 startTimerMemory(mt.makeInterface, mt.memoryInterface);

 DofSet ***allBoundary = new DofSet**[numSub];

 if (domain->solInfo().renum)
   fprintf(stderr," ... Renumbering as Specified       ...\n");

 paralApplyToAll(numSub, localData, &MDAxiData::sendDOFList);

 paralApplyToAll(numSub, localData, &MDAxiData::gatherDOFList, allBoundary);

 stopTimerMemory(mt.makeInterface, mt.memoryInterface);

 startTimerMemory(mt.createDofs, mt.memoryForm);
 int *counterInterf = new int[numSub];
 int **FList = new int*[numSub];
 int **nonZE = new int*[numSub];
 int **subCount = new int*[numSub];
 DComplex ***subPosition = new DComplex**[numSub];
 execParal(numSub, this, &MDAxiDesc::prepareCoarseData, allBoundary,
           counterInterf, FList, nonZE, subCount, subPosition);
 int nP = (numSub>glBCs->numModes) ? numSub : 2*glBCs->numModes+1;
 execParal(nP, this, &MDAxiDesc::makeCoarseData, allBoundary,
           counterInterf, FList, nonZE, subCount, subPosition);
 execParal(numSub, this, &MDAxiDesc::deleteCoarseInfo, FList, nonZE, 
           subCount, subPosition);
 delete[] counterInterf;
 counterInterf = 0;
 stopTimerMemory(mt.createDofs, mt.memoryForm);

 startTimerMemory(mt.makeInterface, mt.memoryInterface);
 delete[] allBoundary;
 allBoundary = 0;
 stopTimerMemory(mt.makeInterface, mt.memoryInterface);

}

void
MDAxiDesc::prepareCoarseData(int iSub, DofSet ***allBoundary,
           int *counterInterf, int **FList, int **nonZE, int **subCount,
           DComplex ***subPosition) {

 localData[iSub]->prepareCoarseData(allBoundary, counterInterf, FList, nonZE,
                  subCount, subPosition);

}


void
MDAxiDesc::deleteCoarseInfo(int iSub, int **FList, int **nonZE, 
           int **subCount, DComplex ***subPosition) {

 localData[iSub]->deleteCoarseInfo(FList, nonZE, subCount, subPosition);

}


void
MDAxiDesc::makeCoarseData(int iP, DofSet ***allBoundary,
           int *counterInterf, int **FList, int **nonZE, int **subCount,
           DComplex ***subPosition) {

 if (numSub>glBCs->numModes) 
   localData[iP]->makeCoarseData(counterInterf, FList, nonZE, subCount,
                  subPosition);
 else
   for (int iSub=0; iSub<numSub; ++iSub)
     localData[iSub]->makeCoarseData(counterInterf, FList, nonZE, subCount,
                      subPosition, iP, iP+1);

}


void
MDAxiDesc::prepareMPCSet(int iP) {

 if (numSub>glBCs->numModes) 
   localData[iP]->prepareMPCSet(); 
 else 
   for (int iSub=0; iSub<numSub; ++iSub)
     localData[iSub]->prepareMPCSet(iP, iP+1);

}


void
MDAxiDesc::getInterfSigns() {

 MatrixTimers &mt = domain->getTimers();
 startTimerMemory(mt.makeInterface, mt.memoryInterface);

 int *wantsToChange = (int *) dbg_alloca(sizeof(int)*numSub);
 int iSub;

 // Self determination of the subdomain sign
 for (iSub=0;iSub<numSub;iSub++) {
    if (iSub % 2  == 0)
      localData[iSub]->InterfSign = 1;
    else
      localData[iSub]->InterfSign = -1;
 }

 int iNeighb;

 int done_flag = 0;

 if (numSub == 1) done_flag = 1;

 while (done_flag == 0) { 
  // Repeat until no change is required

  // Loop over the neighbors to check if a subdomain
  // wants to change its sign
    for (iSub=0;iSub<numSub;iSub++) {

      SComm *sc = localData[iSub]->getSComm();

      int numNeighbors = sc->numNeighb;

      int *InterfaceIndex = localData[iSub]->getInterfaceIndex();

      if(numNeighbors == 0) {
           wantsToChange[iSub] = 1;
           continue;
      }

      wantsToChange[iSub] = 1; // 1: wants to change sign

      for(iNeighb = 0; iNeighb < numNeighbors; iNeighb++) {

        // If iSub and its neighbor do not share faces, skip
        if (InterfaceIndex[iNeighb] == InterfaceIndex[iNeighb+1])
           continue;

        int neighbNum = (sc->subNums)[iNeighb];

        if (localData[iSub]->InterfSign!=localData[neighbNum]->InterfSign) {
          // Found one neighbor with different sign, so NO change sign
          wantsToChange[iSub] = 0;
          break;
        }
      }
    }

    // Finish - unless there are subdomains that want to change sign
    done_flag = 1; 

    for (iSub=0 ; iSub<numSub ; iSub++) 
        if (wantsToChange[iSub]) {
           done_flag = 0;
           break;
        }

    // NOTE: if no one wants to change, then we are done
    if (done_flag) 
       break;

    // Now we will check if a given subdomain has a neighbor that
    // has the same sign and also wants to change sign

    for (iSub=0 ; iSub<numSub ; iSub++) {

      // Only analyze subs who want to change
      if (wantsToChange[iSub] == 0) 
         continue; 

      localData[iSub]->InterfSign *= -1;

      SComm *sc = localData[iSub]->getSComm();

      int numNeighbors = sc->numNeighb;

      int *InterfaceIndex = localData[iSub]->getInterfaceIndex();

      for (iNeighb=0 ; iNeighb<numNeighbors ; iNeighb++) {
         // If iSub and its neighbor do not share faces, skip
         if (InterfaceIndex[iNeighb] == InterfaceIndex[iNeighb+1])
            continue;

         int neighbNum = (sc->subNums)[iNeighb];

         wantsToChange[neighbNum] = 0; 
         // Since I flipped, my neighbors do not need to flip
      }
    }
 }

 stopTimerMemory(mt.makeInterface, mt.memoryInterface);

}


void MDAxiDesc::solve() {

 startTimerMemory(times->preProcess, times->memoryPreProcess);
 preProcess();
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 int MPCsize;
 MPCsize = (glMPCs == 0) ? MPCsize = 0 : 
                           MPCsize = glMPCs->totalNumMPC;

 times->getFetiSolverTime -= getTime();
 solver = new FetiHAxiSolver(numSub, localData, subToSub, MPCsize, 
              &domain->solInfo().getFetiInfo());
 times->getFetiSolverTime += getTime();

 DistrComplexVector rhs(solver->localInfo());
 DistrComplexVector sol(solver->localInfo());

 times->formRhs -= getTime();
 solver->makeStaticLoad(rhs);
 times->formRhs += getTime();

 // Delete elemToNode --- nodeToElem
 paralApplyToAll(numSub, localData, &MDAxiData::deleteConn);

 if (glMPCs != NULL) {
   ComplexVector rhsMPC(MPCsize);
   times->formRhs -= getTime();
   buildRhsMPC(rhsMPC);
   times->formRhs += getTime();
   solver->solveMPC(rhs,rhsMPC,sol);
 }
 else
   solver->solve(rhs,sol);

 startTimerMemory(times->output, times->memoryOutput);
 postProcessing(sol);
 stopTimerMemory(times->output, times->memoryOutput);

 int iSub;
 for (iSub=0; iSub<numSub; ++iSub) {
    times->memoryK += localData[iSub]->getMemoryK();
    times->memoryPrecond += localData[iSub]->getMemoryPrec();
 }

 times->printStaticTimersFetiHAxi(domain->getTimers(), 
        solver->getSolutionTime(), domain->solInfo(), solver->getTimers(),
        geoSource->getCheckFileInfo()[0], domain, glBCs, MPCsize);

}


void
MDAxiDesc::buildRhsMPC(ComplexVector &rhs) {

 int i,j;
 int counter = 0;
 double theta;
 double dot;
 CoordSet &noeuds = domain->getNodes();
 double kappa = localData[0]->getWaveNumber();

 for (i=0; i<glMPCs->numMPCSet; ++i) {
    Node nd = noeuds.getNode(glMPCs->term[i].nnum);
    for (j=0; j<glMPCs->term[i].slices+1; ++j) {
      theta = (j==0) ? glMPCs->term[i].angle1 :
                       glMPCs->term[i].angle1 + (glMPCs->term[i].angle2 -
                       glMPCs->term[i].angle1)*j/glMPCs->term[i].slices;
      dot = nd.x*cos(theta)*glMPCs->term[i].dx + 
            nd.x*sin(theta)*glMPCs->term[i].dy + nd.y*glMPCs->term[i].dz; 
      rhs[counter] = exp(DComplex(0.0,kappa*dot)); 
      rhs[counter++] *= DComplex(glMPCs->term[i].cr,glMPCs->term[i].ci);
    }
 }

}


void
MDAxiDesc::postProcessing(DistrComplexVector &u) {

 fprintf(stderr," Entre dans postprocessing\n");

 double *c = (double *) dbg_alloca(sizeof(double)*glBCs->numSlices+1);
 double *s = (double *) dbg_alloca(sizeof(double)*glBCs->numSlices+1);
 double *coeff = (double *) dbg_alloca(sizeof(double)*glBCs->numModes+1);
 double *scoeff = (double *) dbg_alloca(sizeof(double)*glBCs->numModes+1);

 DComplex result;
 int i, j, k;
 int iInfo, number;
 int maxNode = domain->numNode();

 // Initialize the global solution array
 DComplex **globalsol = new DComplex *[2*glBCs->numModes+1];
 for (i = 0; i<2*glBCs->numModes+1; ++i)
    globalsol[i] = new DComplex[maxNode];
 
 // Merge the solution with the boundary condition
 for (i = 0; i < numSub; ++i) 
    localData[i]->mergeSolution(globalsol, u.subData(i));
 
 for (i=0; i<glBCs->numSlices; ++i) {
    double angle = 2*M_PI*i/glBCs->numSlices;
    c[i] = cos(angle);
    s[i] = sin(angle);
 }

 fprintf(stderr," ... Output Complex Solution        ...\n");

 number = geoSource->getNumOutInfo();
 OutputInfo *sortie = geoSource->getOutputInfo();
 ControlInfo *cinfo = geoSource->getCheckFileInfo();

 for (iInfo=0; iInfo<number; ++iInfo) {

    //int w = sortie[iInfo].width;
    //int p = sortie[iInfo].precision;

    if (sortie[iInfo].type != (OutputInfo::Helmholtz)) {
      continue;   
    }

    if (sortie[iInfo].interval == 1) {
      if ((sortie[iInfo].filptr = fopen(sortie[iInfo].filename,"w")) ==
                 (FILE *) NULL )
         fprintf(stderr,"\n... ERROR: Cannot open %s ...\n", 
                        sortie[iInfo].filename);
      fflush(sortie[iInfo].filptr);
    }

    // Print the real part of the solution

    fprintf(sortie[iInfo].filptr,"Scalar Pressure_Real under load for %s"
            " \n%d\n0.0000000000000\n",cinfo->nodeSetName,
            glBCs->numSlices*maxNode);

    for (i=0; i<glBCs->numSlices; ++i) {
        for (j=0; j<maxNode; ++j) {
            result = globalsol[0][j];
            for (k=1; k<=glBCs->numModes; ++k) {
               int angle = (k*i)%glBCs->numSlices;
               result = result + c[angle]*globalsol[2*k-1][j];
               result = result + s[angle]*globalsol[2*k][j];
            }
//            fprintf(sortie[iInfo].filptr,"% *.*e\n",w,p,real(result) );
            fprintf(sortie[iInfo].filptr,"%1.14e\n",real(result) );
        }
    }

    // Print the imaginary part of the solution

    fprintf(sortie[iInfo].filptr,"Scalar Pressure_Imag under load for %s"
            " \n%d\n0.0000000000000\n",cinfo->nodeSetName,
            glBCs->numSlices*maxNode);

    for (i=0; i<glBCs->numSlices; ++i) {
        for (j=0; j<maxNode; ++j) {
            result = globalsol[0][j];
            for (k=1; k<=glBCs->numModes; ++k) {
               int angle = (k*i)%glBCs->numSlices;
               result = result + c[angle]*globalsol[2*k-1][j];
               result = result + s[angle]*globalsol[2*k][j];
            }
//            fprintf(sortie[iInfo].filptr,"% *.*e\n",w,p,imag(result) );
            fprintf(sortie[iInfo].filptr,"%1.14e\n",imag(result) );
        }
    }

    fflush(sortie[iInfo].filptr);

 }

 for (iInfo=0; iInfo<number; ++iInfo)
    if(sortie[iInfo].interval == 1 && sortie[iInfo].type == 31)
      fclose(sortie[iInfo].filptr);

/*
 // Output the Fourier coefficients of the solution
 //
 // UH - 12-17-99 
 // Put this output file in a more general format

 int len = strlen(cinfo->checkfile);
 char *fileName = (char *) dbg_alloca(sizeof(char)*(len+9));
 strcpy(fileName, cinfo->checkfile);
 strcat(fileName,".fourier");

 FILE *fid = fopen(fileName,"w");

 for (k=0; k<(2*glBCs->numModes+1); ++k) 
   for (j=0; j<maxNode; ++j) {
     fprintf(fid,"%1.7e",real(globalsol[k][j]));
     if (imag(globalsol[k][j]) < 0.0)
        fprintf(fid,"  %1.14e \n",imag(globalsol[k][j]));
     else 
        fprintf(fid,"  %1.14e \n",imag(globalsol[k][j]));
   }
*/

 // UH - 12-17-99
 // 
 // Change this output to see the MAX, MIN, AVERAGE of |u|^2
 //

 for (k=0; k<=glBCs->numModes; ++k) {
     coeff[k] = 0;
     scoeff[k]=0;
     for (j=0; j<maxNode; ++j) {
        if (k==0) {
          coeff[k] += real(globalsol[2*k][j])*real(globalsol[2*k][j]);
          coeff[k] += imag(globalsol[2*k][j])*imag(globalsol[2*k][j]);
        }
        else {
          coeff[k] += real(globalsol[2*k-1][j])*real(globalsol[2*k-1][j]);
          coeff[k] += imag(globalsol[2*k-1][j])*imag(globalsol[2*k-1][j]);
          scoeff[k] += real(globalsol[2*k][j])*real(globalsol[2*k][j]);
          scoeff[k] += imag(globalsol[2*k][j])*imag(globalsol[2*k][j]);
        }
     }
 }

 if (glScatter) 
   buildFFP(globalsol, u);

 verifMPCScat();

 for (i = 0; i<2*glBCs->numModes+1; ++i) {
   delete[] globalsol[i];
   globalsol[i] = 0;
 }
 delete[] globalsol;
 globalsol = 0;

}


void
MDAxiDesc::buildFFP(DComplex **u, DistrComplexVector &sol) {

 if ((glScatter->numFFP<=0) || (glScatter->numElem==0)) {
   fprintf(stderr," ERROR - Inconsistent parameters to build FFP\n");
   return;
 }

 // 3D -> numVector = discretisation de pi en nffp/2

 int nsint = 4*(glScatter->numFFP/4);
 int numPhi = 2*(nsint/4)+1;
 int phiStop = (numPhi+1)/2;
 int i, j;

 DComplex *FFP = new DComplex[numPhi*nsint];
 double (*vectorDir)[3] = new double[phiStop*nsint][3];
 int numSample = 0;

 for (i=0; i<numPhi*nsint; ++i) {
   FFP[i] = DComplex(0.0, 0.0);
 } 

 for (i=0; i<phiStop; ++i) {
   double phi = (numPhi==1) ? 0.0 : M_PI*(-0.5+((double) i)/(numPhi-1.0));
   double c = cos(phi);
   c *= (fabs(c)>1e-15);
   for (j=0; j<nsint; ++j) {
     vectorDir[numSample][0] = c*cos(2*j*M_PI/nsint);
     vectorDir[numSample][1] = c*sin(2*j*M_PI/nsint);
     vectorDir[numSample][2] = sin(phi);
     numSample += 1;
   }
 }

 int nT = threadManager->numThr();
 int nP = (numSub>glBCs->numModes) ? numSub : 2*glBCs->numModes+1;

 int iEle;

 int Fourier;
 int totalFourier = 2*glBCs->numModes+1;        

 // Start the computation of the far field pattern

 if (glScatter->numNode==0) {

   fprintf(stderr," ... FFP on the Obstacle Surface    ...\n");

   // Compute classically the far field pattern on the obstacle

   if (glBCs->numNeu>0) { 

     execParal(nT, this, &MDAxiDesc::ffpAxiNeum, FFP, u, vectorDir); 

   }

   if (glBCs->numDir>0) {

     execParal(nT, this, &MDAxiDesc::ffpAxiDir, FFP, u, vectorDir);

     // Add the contribution of the normal derivative integrated
     DComplex **dudn = new DComplex *[totalFourier];
     for (Fourier=0; Fourier<totalFourier; ++Fourier) {
       dudn[Fourier] = new DComplex[glBCs->numDir]; 
       for (i=0; i<glBCs->numDir; ++i) 
         dudn[Fourier][i] = DComplex(0.0, 0.0);
     }

     execParal(nP, this, &MDAxiDesc::assembleDUDN, &sol, dudn);

     int *pointNumb = (int *) dbg_alloca(sizeof(int)*(glBCs->numDir+1));
     pointNumb[0] = glBCs->numDir;  
     for (i=0; i<glBCs->numDir; ++i)
       pointNumb[i+1] = glBCs->dirBC[i].nnum; 

     execParal(nT, this, &MDAxiDesc::termDUDN, FFP, dudn,vectorDir,pointNumb);

     for (Fourier=0; Fourier<totalFourier; ++Fourier) {
       delete[] dudn[Fourier];
       dudn[Fourier] = 0;
     }
     delete[] dudn;
     dudn = 0;

   }

 }
 else {

   fprintf(stderr," ... FFP on an Inner Surface        ...\n");

   double buffer = -getTime();
   execParal(nT, this, &MDAxiDesc::ffpAxiDir, FFP, u, vectorDir); 
   fprintf(stderr," Fin de la premiere partie : %e\n", 
           (buffer+getTime())/1000.0);

   int nodeNum = domain->numNode();
   int *nodePos = new int[nodeNum];
   int *glToSurf = new int[nodeNum];
   int nodeSurf = 0;
  
   for (i=0; i<nodeNum; ++i) {
     nodePos[i] = 0; 
     glToSurf[i] = -1;
   }

   for (i=0; i<glScatter->numNode; ++i)
     nodePos[glScatter->subNodes[i]] = -1;

   int scatDof = glScatter->scatEle[0]->numDofs();
   int *poinTmp = (int *) dbg_alloca(sizeof(int)*scatDof);

   for (iEle=0; iEle<glScatter->numElem; ++iEle) {
     
     int *dofs = glScatter->scatEle[iEle]->nodes(poinTmp);
     for (i=0; i<scatDof; ++i) {
       if (nodePos[dofs[i]]==0) {
         nodePos[dofs[i]] = 1;
         glToSurf[dofs[i]] = nodeSurf;
         nodeSurf += 1;
       }
     }

   }  

   int *surfToGl = new int[nodeSurf+1];
   surfToGl[0] = nodeSurf;
   for (i=0; i<nodeNum; ++i) {
     if (glToSurf[i]>-1) {
       surfToGl[glToSurf[i]+1] = i;
     }
   }

   // Add the contribution of the normal derivative

   DComplex **dudn = new DComplex *[totalFourier];
   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     dudn[Fourier] = new DComplex[nodeSurf]; 
     for (i=0; i<nodeSurf; ++i) 
       dudn[Fourier][i] = DComplex(0.0, 0.0);
   }

   buffer = -getTime();
   execParal(nP, this, &MDAxiDesc::localDUDN, u, dudn, nodePos, glToSurf);
   fprintf(stderr," Fin de la deuxieme partie : %e\n", 
           (buffer+getTime())/1000.0);

   buffer = -getTime();
   execParal(nT, this, &MDAxiDesc::termDUDN, FFP, dudn, vectorDir, surfToGl); 
   fprintf(stderr," Fin de la troisieme partie : %e\n", 
           (buffer+getTime())/1000.0);

   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     delete[] dudn[Fourier];
     dudn[Fourier] = 0;
   }
   delete[] dudn;
   dudn = 0;

   delete[] nodePos;
   nodePos = 0;
   delete[] glToSurf;
   glToSurf = 0;
   delete[] surfToGl;
   surfToGl = 0;

 }

 // OUTPUT result in FILE
 // 2D -> theta | real part | imaginary part | logarithmic value
 // 3D -> theta | phi | real part | imaginary | logarithmic value

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 int len = strlen(cinfo->checkfile);
 char *fileName = (char *) dbg_alloca(sizeof(char)*(len+8));
 strcpy(fileName, cinfo->checkfile);
 strcat(fileName,".ffp");

 FILE *fffp= fopen(fileName,"w");

 if (numPhi==1) {
   for(i=0;i<nsint;i++) {
     double y = 10.0 * log(2.0*M_PI*abs(FFP[i])*
                abs(FFP[i]))/log(10.0);
     fprintf(fffp,"%e  %.4e  %.4e  %.4e\n", 2*M_PI*i/nsint,
             real(FFP[i]), imag(FFP[i]), y);
   }
 }
 else {
   numSample = 0;
   for (i=0; i<numPhi; ++i) {
     double phi = M_PI*(-0.5+((double) i)/(numPhi-1.0));
     for(j=0; j<nsint; ++j) {
       double y = 10.0 * log(2.0*M_PI*abs(FFP[numSample])*
                  abs(FFP[numSample]))/log(10.0);
       fprintf(fffp,"%e  %e  %.4e  %.4e  %.4e\n", 2*M_PI*j/nsint, phi,
               real(FFP[numSample]),imag(FFP[numSample]),y);
       numSample += 1;
     }
   }
 }

 fclose(fffp);

 // delete FFP & vectorDir
 delete[] FFP;
 FFP = 0; 
 delete[] vectorDir;
 vectorDir = 0;

 fprintf(stderr,"\n --------------------------------------\n");

}


void
MDAxiDesc::termDUDN(int iT, DComplex *FFP, DComplex* *dudn,  
                    double (*vectorDir)[3], int *pointNumb) {

 double kappa = localData[0]->getWaveNumber();
 int totalFourier = 2*glBCs->numModes+1;
 int i, j, k;
 int Fourier;
 DComplex coeff = DComplex(0.25/M_PI, 0.0);

 CoordSet &noeuds = domain->getNodes();

 int nT = threadManager->numThr();
 int nsint = 4*(glScatter->numFFP/4);
 int numPhi = 2*(nsint/4)+1;
 int phiStop = (numPhi+1)/2;

// No need for mutual exclusion nor atomic operations

 for (k=iT; k<phiStop; k+=nT) {
   double phi = (numPhi==1) ? 0.0 : M_PI*(-0.5+((double) k)/(numPhi-1.0));
   double rho = cos(phi);
   rho *= (rho>1e-15);
   for (i=0; i<pointNumb[0]; ++i) {
     Node nd = noeuds.getNode(pointNumb[i+1]);
     double eta = kappa*nd.y*vectorDir[k*nsint][2];
     DComplex C = DComplex(cos(eta),-sin(eta)); 
     double rhoX = fabs(nd.x)*rho;
     for (j=0; j<nsint; ++j) {
       int position = k*nsint+j;
       double ratio = (rhoX<1e-15) ? 1.0 : vectorDir[position][0]*nd.x/rhoX;
       ratio = (ratio>1.0)  ?  1.0 : ratio;
       ratio = (ratio<-1.0) ? -1.0 : ratio;
       double alpha = acos(ratio);
       if (vectorDir[position][1]<0.0)
         alpha *= -1.0;
       for (Fourier=0; Fourier<totalFourier; ++Fourier) {
         DComplex expDir = coefExpDir(Fourier, kappa, -rhoX, alpha)*C;
         if ((real(expDir)==0.0) && (imag(expDir)==0.0)) 
           continue;
         FFP[position] += expDir*dudn[Fourier][i]*coeff;
         if (k<phiStop-1) {
           int position2 = ((j+nsint/2)%nsint) + nsint*(numPhi-1-k); 
           FFP[position2] += DComplex(real(expDir),-imag(expDir))*
                             dudn[Fourier][i]*coeff;
         } 
       }
     } 
   } 
 }

}


void
MDAxiDesc::ffpAxiNeum(int iT, DComplex *FFP, DComplex **u, 
                      double (*vectorDir)[3]) {

 double kappa = localData[0]->getWaveNumber();
 int totalFourier = 2*glBCs->numModes+1;

 double *waveDirection = (double *) dbg_alloca(sizeof(double)*3);
 waveDirection[0] = glBCs->dirX;
 waveDirection[1] = glBCs->dirY;
 waveDirection[2] = glBCs->dirZ;

 int nT = threadManager->numThr();

 int Fourier;
 int iEle;
 int i, j, k;

 CoordSet &noeuds = domain->getNodes();

 int scatDof = glScatter->scatEle[0]->numDofs();
 int *poinTmp = (int*) dbg_alloca(sizeof(int)*scatDof);

 DComplex **uel = (DComplex **) dbg_alloca(sizeof(DComplex*)*totalFourier);
 for (Fourier=0; Fourier<totalFourier; ++Fourier) 
   uel[Fourier] = (DComplex *) dbg_alloca(sizeof(DComplex)*scatDof);

 int nsint = 4*(glScatter->numFFP/4);
 int numPhi = 2*(nsint/4)+1;
 int phiStop = (numPhi+1)/2;

 DComplex *bufferFFP = (DComplex *) dbg_alloca(sizeof(DComplex)*2*nsint);

 for (iEle=iT; iEle<glScatter->numElem; iEle+=nT) {

   LineAxiSommer *ele1=dynamic_cast<LineAxiSommer *>(glScatter->scatEle[iEle]);
   Line2AxiSommer *ele2 = dynamic_cast<Line2AxiSommer *>
                          (glScatter->scatEle[iEle]);

   if ((ele1==0) && (ele2==0)) {
     fprintf(stderr,"Element chosen to build FFP non axisymmetric."
                    " Exiting\n");
     return;
   }

   int *dofs = (ele1) ? ele1->nodes(poinTmp) : ele2->nodes(poinTmp);

   for (Fourier=0; Fourier<totalFourier; ++Fourier) 
     for (i=0; i<scatDof; ++i) 
           uel[Fourier][i] = u[Fourier][dofs[i]];

   for (k=0; k<phiStop; ++k) {

     int kPhi = (k+iT)%phiStop;

     for (j=0; j<nsint; ++j) {
       bufferFFP[2*j] = DComplex(0.0, 0.0);
       bufferFFP[2*j+1] = DComplex(0.0, 0.0);
     }

     if (ele1) { 
       ele1->ffpAxiNeum(nsint, bufferFFP, noeuds, uel, kappa,
             vectorDir+kPhi*nsint, waveDirection, totalFourier);
     }
     else {
       ele2->ffpAxiNeum(nsint, bufferFFP, noeuds, uel, kappa,
             vectorDir+kPhi*nsint, waveDirection, totalFourier);
     }

#if defined(sgi) && ! defined(_OPENMP)
     ussetlock(allocLock);
     for (j=0; j<nsint; ++j) {
       FFP[kPhi*nsint+j] += bufferFFP[2*j];
       if (kPhi<phiStop-1) {
         int position2 = ((j+nsint/2)%nsint) + nsint*(numPhi-1-kPhi); 
         FFP[position2] += bufferFFP[2*j+1];
       } 
     } 
     usunsetlock(allocLock);
#else
     for (j=0; j<nsint; ++j) {
       FFP[kPhi*nsint+j] += bufferFFP[2*j];
       if (kPhi<phiStop-1) {
         int position2 = ((j+nsint/2)%nsint) + nsint*(numPhi-1-kPhi); 
         FFP[position2] += bufferFFP[2*j+1];
       } 
     } 
#endif

   }
 }

}


void
MDAxiDesc::ffpAxiDir(int iT, DComplex *FFP, DComplex **u, 
                      double (*vectorDir)[3]) {

 double kappa = localData[0]->getWaveNumber();
 int totalFourier = 2*glBCs->numModes+1;
 int nT = threadManager->numThr();

 int Fourier;
 int iEle;
 int i, j, k;

 CoordSet &noeuds = domain->getNodes();

 int scatDof = glScatter->scatEle[0]->numDofs();
 int *poinTmp = (int*) dbg_alloca(sizeof(int)*scatDof);

 DComplex **uel = (DComplex **) dbg_alloca(sizeof(DComplex*)*totalFourier);
 for (Fourier=0; Fourier<totalFourier; ++Fourier) 
   uel[Fourier] = (DComplex *) dbg_alloca(sizeof(DComplex)*scatDof);

 int nsint = 4*(glScatter->numFFP/4);
 int numPhi = 2*(nsint/4)+1;
 int phiStop = (numPhi+1)/2;

 DComplex *bufferFFP = (DComplex *) dbg_alloca(sizeof(DComplex)*2*nsint);

 double *waveDirection = (double *) dbg_alloca(sizeof(double)*3);
 if (glBCs->numDir>0) {
   waveDirection[0] = glBCs->dirX;
   waveDirection[1] = glBCs->dirY;
   waveDirection[2] = glBCs->dirZ;
 }
 else {
   waveDirection[0] = glBCs->dirX;
   waveDirection[1] = glBCs->dirY;
   waveDirection[2] = glBCs->dirZ;
 }

 for (iEle=iT; iEle<glScatter->numElem; iEle+=nT) {

   LineAxiSommer *ele1=dynamic_cast<LineAxiSommer *>(glScatter->scatEle[iEle]);
   Line2AxiSommer *ele2 = dynamic_cast<Line2AxiSommer *>
                          (glScatter->scatEle[iEle]);
 
   if ((ele1==0) && (ele2==0)) {
     fprintf(stderr,"Element chosen to build FFP non axisymmetric."
                        " Exiting\n");
     return;
   }

   int *dofs = (ele1) ? ele1->nodes(poinTmp) : ele2->nodes(poinTmp);

   for (Fourier=0; Fourier<totalFourier; ++Fourier) 
     for (i=0; i<scatDof; ++i) {
       uel[Fourier][i] = u[Fourier][dofs[i]];
     }

   for (k=0; k<phiStop; ++k) {

     int kPhi = (k+iT)%phiStop;

     for (j=0; j<nsint; ++j) {
       bufferFFP[2*j] = DComplex(0.0, 0.0);
       bufferFFP[2*j+1] = DComplex(0.0, 0.0);
     }
 
     if (ele1) {
       ele1->ffpAxiDir(nsint, bufferFFP, noeuds, uel, kappa,
            vectorDir+kPhi*nsint, waveDirection, totalFourier);
     }
     else { 
       ele2->ffpAxiDir(nsint, bufferFFP, noeuds, uel, kappa,
            vectorDir+kPhi*nsint, waveDirection, totalFourier);
     }

#if defined(sgi) && ! defined(_OPENMP)
     ussetlock(allocLock);
     for (j=0; j<nsint; ++j) {
       FFP[kPhi*nsint+j] += bufferFFP[2*j];
       if (kPhi<phiStop-1) {
         int position2 = ((j+nsint/2)%nsint) + nsint*(numPhi-1-kPhi); 
         FFP[position2] += bufferFFP[2*j+1];
       } 
     } 
     usunsetlock(allocLock);
#else 
     for (j=0; j<nsint; ++j) {
       FFP[kPhi*nsint+j] += bufferFFP[2*j];
       if (kPhi<phiStop-1) {
         int position2 = ((j+nsint/2)%nsint) + nsint*(numPhi-1-kPhi);
         FFP[position2] += bufferFFP[2*j+1];
       }
     }
#endif

   }
 }

}


void
MDAxiDesc::assembleDUDN(int iP, DistrComplexVector *sol, DComplex **dudn) {

 if (numSub>glBCs->numModes) 
   localData[iP]->assembleDUDN(sol->subData(iP), dudn);
 else
   for (int iSub=0; iSub<numSub; ++iSub) {
     localData[iSub]->assembleDUDN(sol->subData(iSub), dudn, iP, iP+1);
   }

}


void
MDAxiDesc::localDUDN(int iP, DComplex **u, DComplex **dudn, int *nodePos, 
           int *glToSurf) {

 if (numSub>glBCs->numModes) 
   localData[iP]->localDUDN(u, dudn, nodePos, glToSurf);
 else
   for (int iSub=0; iSub<numSub; ++iSub) {
     localData[iSub]->localDUDN(u, dudn, nodePos, glToSurf, iP, iP+1);
   }

}


void
MDAxiDesc::verifMPCScat() {

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 int len = strlen(cinfo->checkfile);
 char *fileName = (char *) dbg_alloca(sizeof(char)*(len+8));
 strcpy(fileName, cinfo->checkfile);
 strcat(fileName,".verif");

 int i;

 FILE *fid = fopen(fileName,"w");

 if ((glScatter) && (glScatter->numElem>0)) {

   fprintf(fid,"Elements ElScatterer using %s\n", cinfo->nodeSetName);
   int scatDof = glScatter->scatEle[0]->numDofs();
   int *pointer = (int *) dbg_alloca(sizeof(int)*scatDof);

   for (i=0; i<glScatter->numElem; ++i) {
     int *nn = glScatter->scatEle[i]->nodes(pointer); 
     fprintf(fid," %d 4 %d %d %d \n", i+1, nn[0]+1, nn[0]+1, nn[1]+1); 
   }

 }
 
 int maxNode = domain->numNode();

 if ((glScatter) && (glScatter->numNode>0)) {
    
   fprintf(fid,"Scalar NodeEliminated under CHECKING for %s\n",
           cinfo->nodeSetName);     
   fprintf(fid,"%d\n0.0000000000000\n", glBCs->numSlices*maxNode);

   double *tmp = new double[glBCs->numSlices*maxNode];
   for (i=0; i<glBCs->numSlices*maxNode; ++i)
     tmp[i] = 0.0;
   for (i=0; i<glScatter->numNode; ++i)
     tmp[glScatter->subNodes[i]] = -1.0;
   for (i=0; i<glBCs->numSlices*maxNode; ++i)
     fprintf(fid,"%e\n", tmp[i]);

   delete[] tmp;
   tmp = 0;

 }

 if (glMPCs) {

   fprintf(fid,"Scalar NodeMPC under CHECKING for %s\n",
           cinfo->nodeSetName);
   fprintf(fid,"%d\n0.0000000000000\n", glBCs->numSlices*maxNode);

   double *tmp = new double[glBCs->numSlices*maxNode];
   for (i=0; i<glBCs->numSlices*maxNode; ++i)
     tmp[i] = -1.0;
   for (i=0; i<glMPCs->numMPCSet; ++i)
     tmp[glMPCs->term[i].nnum] = glMPCs->term[i].slices+1;
   for (i=0; i<glBCs->numSlices*maxNode; ++i)
     fprintf(fid,"%e\n", tmp[i]);

   delete[] tmp;
   tmp = 0;

 } 

 fclose(fid);

}

